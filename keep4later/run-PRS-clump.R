
bigassertr::assert_dir("results-clump")
bigassertr::assert_dir("log")

library(dplyr)
files <- tibble(pheno = c(paste0("t1d_", c("affy", "illu")),
                          paste0("brca_", c("onco", "icogs")),
                          "cad", "prca", "mdd")) %>%
  tidyr::expand_grid(
    qc = c("qc1", "qc2"),
    name_corr = c("data/corr-large/clump_adj_chr",
                  "data/corr-large/clump_adj_with_blocks_chr")) %>%
  print()

library(future.batchtools)
NCORES <- 30
plan(batchtools_slurm(resources = list(
  t = "24:00:00", c = NCORES + 2, mem = "150g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pwalk(files, function(pheno, qc, name_corr) {

  # pheno <- "t1d_illu"
  # qc <- "qc1"
  # name_corr <- "data/corr-large/clump_adj_with_blocks_chr"

  gwas_file <- paste0("data/sumstats-clump/", pheno, "_", qc, ".rds")
  sumstats <- readRDS(gwas_file)
  basename <- file.path("results-clump", sub("\\.rds$", "", basename(gwas_file)))
  if (grepl("blocks", name_corr, fixed = TRUE))
    basename <- paste0(basename, "_adj_with_blocks")

  library(bigsnpr)
  CHR <- as.integer(readRDS("tmp-data/fake_for_clump.rds")$map$chromosome)

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
    file = paste0(basename, "_LDpred2", ".rds"))

  (h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4))
  (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(TRUE, FALSE)))

  runonce::save_run(
    snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES),
    file = paste0(basename, "_LDpred2-low-h2", ".rds"))


  # lassosum2
  runonce::save_run(
    snp_lassosum2(corr, df_beta, ncores = NCORES),
    file = paste0(basename, "_lassosum2", ".rds"))


  # LDpred2-auto
  runonce::save_run(
    snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, ncores = NCORES,
                     vec_p_init = seq_log(1e-4, 0.5, 30),
                     allow_jump_sign = TRUE, shrink_corr = 1),
    file = paste0(basename, "_LDpred2-auto", ".rds"))

  runonce::save_run(
    snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, ncores = NCORES,
                     vec_p_init = seq_log(1e-4, 0.5, 30),
                     allow_jump_sign = FALSE, shrink_corr = 0.95),
    file = paste0(basename, "_LDpred2-auto-rob", ".rds"))
})
