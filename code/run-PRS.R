library(dplyr)
grid <- tibble(pheno = c(paste0("t1d_", c("affy", "illu")),
                         paste0("brca_", c("onco", "icogs")),
                         "cad", "prca", "mdd")) %>%
  mutate(whichN = NA) %>%
  add_row(pheno = "vitaminD", whichN = c("trueN", "maxN")) %>%
  tidyr::expand_grid(
    qc = c("noqc", "qc1", "qc2"),
    name_corr = c("data/corr/chr", "data/corr/adj_with_blocks_chr")) %>%
  print()

library(future.batchtools)
NCORES <- 30
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "120g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results")
bigassertr::assert_dir("results/sumstats")

furrr::future_pwalk(grid, function(pheno, whichN, qc, name_corr) {

  # pheno <- "vitaminD"
  # whichN <- "trueN"
  # qc <- "qc1"
  # name_corr <- "data/corr/adj_with_blocks_chr"

  gwas_file <- paste0("data/sumstats/", pheno, "_", qc, ".rds")
  sumstats <- readRDS(gwas_file)
  if (identical(whichN, "maxN"))
    sumstats$n_eff <- max(sumstats$n_eff)

  suffix_N <- `if`(is.na(whichN), "", paste0("_", whichN))
  basename <- file.path("results/sumstats", sub("\\.rds$", suffix_N, basename(gwas_file)))
  if (grepl("blocks", name_corr, fixed = TRUE))
    basename <- paste0(basename, "_adj_with_blocks")

  if (file.exists(paste0(basename, "_LDpred2-auto-rob.rds")))
    return(NULL)

  library(bigsnpr)
  CHR <- snp_attach("data/UKBB_HM3_val.rds")$map$chromosome

  tmp <- tempfile(tmpdir = "tmp-data/sfbm/")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

  for (chr in 1:22) {

    print(chr)

    ## indices in 'sumstats'
    ind.chr <- which(sumstats$chr == chr)
    ## indices in 'G'
    ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
    ## indices in 'corr'
    ind.chr3 <- match(ind.chr2, which(CHR == chr))

    corr0 <- readRDS(paste0(name_corr, chr, ".rds"))
    ld_chr <- Matrix::colSums(corr0^2)
    df_beta_chr <- sumstats[ind.chr, c("beta", "beta_se", "n_eff")]

    if (chr == 1) {
      df_beta <- df_beta_chr
      ld <- ld_chr[ind.chr3]
      corr <- as_SFBM(corr0[ind.chr3, ind.chr3], tmp, compact = TRUE)
    } else {
      df_beta <- rbind(df_beta, df_beta_chr)
      ld <- c(ld, ld_chr[ind.chr3])
      corr$add_columns(corr0[ind.chr3, ind.chr3], nrow(corr))
    }
  }


  # Run LD score regression
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))
  h2_est <- ldsc[["h2"]]


  # LDpred2-grid
  (h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
  (p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
  (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(TRUE, FALSE)))

  runonce::save_run(
    snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES),
    file = paste0(basename, "_LDpred2.rds"))

  (h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4))
  (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(TRUE, FALSE)))

  runonce::save_run(
    snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES),
    file = paste0(basename, "_LDpred2-low-h2.rds"))


  # lassosum2
  runonce::save_run(
    snp_lassosum2(corr, df_beta, ncores = NCORES),
    file = paste0(basename, "_lassosum2.rds"))


  # LDpred2-auto
  runonce::save_run(
    snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, ncores = NCORES,
                     vec_p_init = seq_log(1e-4, 0.5, 30),
                     allow_jump_sign = TRUE, shrink_corr = 1),
    file = paste0(basename, "_LDpred2-auto.rds"))

  runonce::save_run(
    snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, ncores = NCORES,
                     vec_p_init = seq_log(1e-4, 0.5, 30),
                     allow_jump_sign = FALSE, shrink_corr = 0.95),
    file = paste0(basename, "_LDpred2-auto-rob.rds"))
})


#### Run PRS-CS ####

library(dplyr)
grid <- tibble(pheno = c(paste0("t1d_", c("affy", "illu")),
                         paste0("brca_", c("onco", "icogs")),
                         "cad", "prca", "mdd", "vitaminD")) %>%
  tidyr::expand_grid(qc = c("noqc", "qc1", "qc2"), chr = 1:22) %>%
  filter(!file.exists(paste0("results_prscs/sumstats/", pheno, "_", qc, "_chr", chr, ".rds"))) %>%
  print()

library(future.batchtools)
NCORES <- 10
plan(batchtools_slurm(workers = 1000, resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "50g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results_prscs")
bigassertr::assert_dir("results_prscs/sumstats")

furrr::future_pwalk(grid, function(pheno, qc, chr) {

  # pheno <- "vitaminD"
  # qc <- "qc1"
  # chr <- 22

  # prepare sumstats file
  sumstats <- readRDS(paste0("data/sumstats/", pheno, "_", qc, ".rds"))
  sumstats_file <- tempfile(tmpdir = "tmp-data", fileext = ".txt")
  on.exit(file.remove(sumstats_file), add = TRUE)
  sumstats %>%
    transmute(SNP = paste0("snp", `_NUM_ID_`), A1 = "A", A2 = "C", BETA = beta,
              P = pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE)) %>%
    bigreadr::fwrite2(sumstats_file, sep = "\t") %>%
    readLines(n = 5) %>%
    writeLines()

  # run PRS-CS for one chromosome
  dir <- "data/corr/ukbb_for_prscs"
  prefix <- tempfile(tmpdir = "tmp-data")
  res_file <- paste0(prefix, "_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt")
  on.exit(file.remove(res_file), add = TRUE)

  system(glue::glue(
    "OMP_NUM_THREADS=", NCORES,
    " python3 PRScs/PRScs.py",
    " --ref_dir={dir}",
    " --bim_prefix={dir}/for_prscs",
    " --sst_file={sumstats_file}",
    " --n_gwas={as.integer(max(sumstats$n_eff))}",
    " --chrom={chr}",
    " --out_dir={prefix}"
  ))

  saveRDS(bigreadr::fread2(res_file),
          paste0("results_prscs/sumstats/", pheno, "_", qc, "_chr", chr, ".rds"))
})
