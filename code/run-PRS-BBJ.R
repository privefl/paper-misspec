ALL_CORR <- c(
  corr_JPN               = "data/corr_JPN/chr",
  corr_JPN_with_blocks   = "data/corr_JPN/adj_with_blocks_chr",
  corr_1000G             = "data/corr_1000G_EAS/chr",
  corr_1000G_with_blocks = "data/corr_1000G_EAS/adj_with_blocks_chr",
  corr_UKBB              = "data/corr_UKBB_EAS/chr",
  corr_UKBB_with_blocks  = "data/corr_UKBB_EAS/adj_with_blocks_chr"
)

library(dplyr)
grid <- tibble(pheno = c("height", "systolic_bp", "hdl_cholesterol", "bmi")) %>%
  tidyr::expand_grid(name_corr = names(ALL_CORR)) %>%
  filter(!file.exists(paste0("results/sumstats_BBJ/", pheno, "_", name_corr,
                             "_LDpred2-auto-rob.rds"))) %>%
  print()

library(future.batchtools)
NCORES <- 14
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "60g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results/sumstats_BBJ")

furrr::future_pwalk(grid, function(pheno, qc, name_corr) {

  library(bigsnpr)
  corr_path <- ALL_CORR[[name_corr]]
  grp <- readRDS(file.path(dirname(corr_path), "all_final_grp.rds"))
  if (is.null(grp$ind)) {
    CHR <- snp_attach("data/UKBB_HM3_val.rds")$map$chromosome
    ind_keep <- split(seq_along(CHR), CHR)
  } else {
    ind_keep <- grp$ind
  }

  sumstats <- readRDS(paste0("data/sumstats_BBJ/", pheno, ".rds")) %>%
    filter(`_NUM_ID_` %in% unlist(ind_keep))
  basename <- paste0("results/sumstats_BBJ/", pheno, "_", name_corr)

  tmp <- tempfile(tmpdir = "tmp-data/sfbm/")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

  for (chr in 1:22) {

    print(chr)

    ## indices in 'sumstats'
    ind.chr <- which(sumstats$chr == chr)
    ## indices in 'G'
    ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
    ## indices in 'corr'
    ind.chr3 <- match(ind.chr2, ind_keep[[chr]])

    corr0 <- readRDS(paste0(corr_path, chr, ".rds"))
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
ALL_CORR <- c(
  corr_JPN   = "data/corr_JPN/ukbb_for_prscs/ldblk_ukbb_chr1.hdf5",
  corr_1000G = "data/corr_1000G_EAS/1kg_for_prscs/ldblk_1kg_chr1.hdf5",
  corr_UKBB  = "data/corr_UKBB_EAS/ukbb_for_prscs/ldblk_ukbb_chr1.hdf5"
)
grid <- tibble(pheno = c("height", "systolic_bp", "hdl_cholesterol", "bmi")) %>%
  tidyr::expand_grid(name_corr = names(ALL_CORR), chr = 1:22) %>%
  filter(!file.exists(paste0("results_prscs/sumstats_BBJ/", pheno, "_", name_corr, "_chr", chr, ".rds"))) %>%
  print()

library(future.batchtools)
NCORES <- 10
plan(batchtools_slurm(workers = 1000, resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "50g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results_prscs/sumstats_BBJ")

furrr::future_pwalk(grid, function(pheno, name_corr, chr) {

  sumstats <- readRDS(paste0("data/sumstats_BBJ/", pheno, ".rds"))
  sumstats_file <- tempfile(tmpdir = "tmp-data", fileext = ".txt")
  on.exit(file.remove(sumstats_file), add = TRUE)
  sumstats %>%
    transmute(SNP = paste0("snp", `_NUM_ID_`), A1 = "A", A2 = "C", BETA = beta,
              P = pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE)) %>%
    bigreadr::fwrite2(sumstats_file, sep = "\t") %>%
    readLines(n = 5) %>%
    writeLines()

  # run PRS-CS for one chromosome
  dir <- dirname(ALL_CORR[name_corr])
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
          paste0("results_prscs/sumstats_BBJ/", pheno, "_", name_corr, "_chr", chr, ".rds"))
})
