
#### Run LDpred2 ####

run_ldpred2 <- function(df_beta, corr) {

  # LDSc reg
  corr0 <- readRDS("tmp-data/corr0_simu_val.rds")
  print(ldsc <- snp_ldsc2(corr0, df_beta))
  h2_est <- ldsc[["h2"]]

  # LDpred-inf
  beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
  pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test)
  print(r21 <- cor(pred_inf, y[ind.test])**2)

  ## LDpred2-grid
  (h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
  (p_seq <- signif(seq_log(1e-4, 1, length.out = 15), 2))
  (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = FALSE))

  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
  pred_grid <- big_prodMat(G, beta_grid, ncores = NCORES)
  score.val <- apply(pred_grid[ind.val, ], 2, cor, y = y[ind.val])
  pred_gibbs <- pred_grid[ind.test, which.max(score.val)]
  print(r22 <- cor(pred_gibbs, y[ind.test])**2)

  # LDpred2-grid (including small values for h2)
  (h2_seq2 <- round(h2_est * c(0.01, 0.1, 0.3, 0.7, 1, 1.4), 4))
  (params2 <- expand.grid(p = p_seq, h2 = h2_seq2, sparse = FALSE))

  beta_grid2 <- snp_ldpred2_grid(corr, df_beta, params2, ncores = NCORES)
  pred_grid2 <- big_prodMat(G, beta_grid2, ncores = NCORES)
  score.val2 <- apply(pred_grid2[ind.val, ], 2, cor, y = y[ind.val])
  pred_gibbs2 <- pred_grid2[ind.test, which.max(score.val2)]
  print(r23 <- cor(pred_gibbs2, y[ind.test])**2)

  # LDpred-auto
  auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                           shrink_corr = 1, allow_jump_sign = TRUE)
  print(c(auto[[1]]$h2_est, auto[[1]]$p_est))
  pred_auto <- big_prodVec(G, auto[[1]]$beta_est, ind.row = ind.test)
  print(r24 <- cor(pred_auto, y[ind.test])**2)

  # LDpred-auto (with new robustness changes)
  auto2 <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                            shrink_corr = 0.9, allow_jump_sign = FALSE)
  print(c(auto2[[1]]$h2_est, auto2[[1]]$p_est))
  pred_auto2 <- big_prodVec(G, auto2[[1]]$beta_est, ind.row = ind.test)
  print(r25 <- cor(pred_auto2, y[ind.test])**2)

  c(r21, r22, r23, r24, r25)
}

# ldpred2 <- run_ldpred2(df_beta)
# mis_ldpred2 <- run_ldpred2(df_beta2)


#### Run lassosum2 ####

run_lassosum2 <- function(df_beta, corr) {
  beta_lassosum <- snp_lassosum2(corr, df_beta, ncores = NCORES)
  pred_grid <- big_prodMat(G, beta_lassosum, ncores = NCORES)
  score.val <- apply(pred_grid[ind.val, ], 2, cor, y = y[ind.val])
  pred_lassosum <- pred_grid[ind.test, which.max(score.val)]
  print(cor(pred_lassosum, y[ind.test])**2)
}

# lassosum2 <- run_lassosum2(df_beta)
# mis_lassosum2 <- run_lassosum2(df_beta2)


#### Run C+T ####

run_CT <- function(df_beta) {

  CHR <- as.integer(data$map$chromosome)
  POS <- data$map$physical.pos
  all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.val,
                                grid.base.size = 200,
                                grid.thr.r2 = c(0.05, 0.2, 0.8),
                                lpS = df_beta$lpval, ncores = NCORES)
  attr(all_keep, "grid")

  multi_PRS <- snp_grid_PRS(G, all_keep, df_beta$beta, df_beta$lpval,
                            ind.row = ind.val, n_thr_lpS = 50, ncores = NCORES)
  dim(multi_PRS)  ## 150 C+T scores

  library(dplyr)
  grid2 <- attr(all_keep, "grid") %>%
    mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = row_number()) %>%
    tidyr::unnest(cols = "thr.lp")
  s <- nrow(grid2)
  grid2$score <- big_univLinReg(multi_PRS, y[ind.val])$score
  max_prs <- grid2 %>% arrange(desc(score)) %>% slice(1:10) %>% print() %>% slice(1)

  ind.keep <- unlist(purrr::map(all_keep, max_prs$id))
  pred_CT <- snp_PRS(G, df_beta$beta[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
                     lpS.keep = df_beta$lpval[ind.keep], thr.list = max_prs$thr.lp)
  print(cor(pred_CT, y[ind.test])[[1]]**2)
}

# CT <- run_CT(df_beta)


#### Run lassosum ####

run_lassosum <- function(df_beta) {

  t <- with(df_beta, beta / beta_se)
  n <- df_beta$n_eff

  library(lassosum)
  doParallel::registerDoParallel(cl <- parallel::makeCluster(NCORES))
  system.time(
    out <- lassosum.pipeline(
      cor = t / sqrt(n - 2 + t^2),
      snp = data$map$marker.ID,
      A1 = data$map$allele1,
      A2 = data$map$allele2,
      exclude.ambiguous = FALSE,
      test.bfile = "tmp-data/simu_chr22",
      LDblocks = "EUR.hg19",
      cluster = cl,
      destandardize = TRUE
    )
  ) # < 3 min
  parallel::stopCluster(cl)

  all_score <- apply(do.call("cbind", out$pgs), 2, cor, y = y[ind.val])
  best_beta <- do.call("cbind", out$beta)[, which.max(all_score)]
  pred_lassosum <- big_prodVec(G, best_beta, ind.row = ind.test)
  print(cor(pred_lassosum, y[ind.test])**2)
}

# lassosum <- run_lassosum(df_beta)
# mis_lassosum <- run_lassosum(df_beta2)


#### Run SBayesR -> do not run because always diverged ####

# gctb <- "../paper-ldpred2/tmp-data/gctb_2.02_Linux/gctb"
#
# # Compute LD
# if (!file.exists(paste0("tmp-data/ldm_22.ldm.shrunk.bin"))) {
#   system.time(
#     system(glue::glue(
#       "{gctb} --bfile tmp-data/simu_chr22",
#       " --make-shrunk-ldm",
#       " --out tmp-data/ldm_22"
#     ))
#   ) # 18 min
# }
#
# # Compute SBayesR
# obj.bed <- bigsnpr::bed("tmp-data/simu_chr22.bed")
# af <- bigsnpr::bed_MAF(obj.bed, ncores = NCORES)$af
#
# tmp <- tempfile(tmpdir = "tmp-data", fileext = ".ma")
# library(dplyr)
# df_beta %>%
#   bind_cols(data$map) %>%
#   transmute(SNP = rsid, A1 = allele1, A2 = allele2, freq = af,
#             b = beta, se = beta_se, p = 10^-lpval, N = n_eff) %>%
#   bigreadr::fwrite2(tmp, sep = " ") %>%
#   readLines(n = 5) %>%
#   writeLines()
#
# res_file <- "tmp-data/sbayesr_chr22"
#
# system(glue::glue(
#   gctb,
#   " --ldm tmp-data/ldm_{chr}.ldm.shrunk",
#   " --sbayes R --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1",
#   # " --sbayes R --pi 0.9,0.1 --gamma 0.0,0.1",
#   # " --p-value 0.4 --rsq 0.95",
#   " --gwas-summary {tmp}",
#   " --chain-length 10000 --burn-in 2000",
#   " --out {res_file} --out-freq 100"
# ))
#
# file.remove(tmp)
#
# library(dplyr)
# head(res_sbayesr <- bigreadr::fread2(paste0(res_file, ".snpRes")))
#
# ind <- match(res_sbayesr$Name, data$map$rsid)
# stopifnot(!anyNA(ind))
# stopifnot(all.equal(data$map$allele1[ind], res_sbayesr$A1))
# stopifnot(all.equal(data$map$allele2[ind], res_sbayesr$A2))
#
# pred_sbayesr <- big_prodVec(G, res_sbayesr$A1Effect, ind.row = ind.test,
#                             ind.col = ind)
# cor(pred_sbayesr, y[ind.test])**2


#### Run PRS-CS -> do not run because not enough variants matched with their LD ref ####

# run_prscs <- function(df_beta) {
#
#   tmp <- tempfile(tmpdir = "tmp-data", fileext = ".txt")
#   df_beta %>%
#     bind_cols(data$map) %>%
#     transmute(SNP = rsid, A1 = allele1, A2 = allele2,
#               BETA = beta, P = 10^-lpval) %>%
#     bigreadr::fwrite2(tmp, sep = "\t") %>%
#     readLines(n = 5) %>%
#     writeLines()
#
#   on.exit(file.remove(tmp), add = TRUE)
#
#   prefix <- paste0("tmp-data/TMP_PRSCS", ic)
#
#   PHI <- c(NA, 10^(-4:0))
#   prscs <- "PRScs/PRScs.py"
#   for (phi in PHI) {
#     system(glue::glue(
#       "OMP_NUM_THREADS=", NCORES[[1]],
#       " python3 {prscs}",
#       " --ref_dir=ldblk_1kg_eur",
#       " --bim_prefix=tmp-data/simu_chr22",
#       " --sst_file={tmp}",
#       " --n_gwas={max(df_beta$n_eff)}",
#       if (is.na(phi)) "" else " --phi={phi}",
#       " --chrom=22",
#       " --out_dir={prefix}"
#     ))
#   }
#   on.exit(unlink(paste0(prefix, "*")), add = TRUE)
#
#   get_pred <- function(phi) {
#     file <- paste0(prefix, "_pst_eff_a1_b0.5_phi", phi, "_chr22.txt")
#     res <- bigreadr::fread2(file)
#     betas <- rep(0, ncol(G))
#     betas[match(res$V2, data$map$rsid)] <- res$V6
#     big_prodVec(G, betas, ncores = NCORES)
#   }
#
#   all_pred <- sapply(sprintf("%.0e", PHI[-1]), get_pred)
#   ind.best <- which.max(apply(all_pred[ind.val, ], 2, cor, y = y[ind.val]))
#   pred_prscs <- all_pred[ind.test, ind.best]
#   print(r21 <- cor(pred_prscs, y[ind.test])**2)
#
#   pred_prscs_auto <- get_pred("auto")[ind.test]
#   print(r22 <- cor(pred_prscs_auto, y[ind.test])**2)
#
#   c(r21, r22)
# }

# prscs <- run_prscs(df_beta)
# mis_prscs <- run_prscs(df_beta2)
