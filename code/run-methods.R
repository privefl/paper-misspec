
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
  print(score.val <- apply(pred_grid[ind.val, ], 2, cor, y = y[ind.val]))
  pred_lassosum <- pred_grid[ind.test, which.max(score.val)]
  print(cor(pred_lassosum, y[ind.test])**2)
}

# lassosum2 <- run_lassosum2(df_beta)
# mis_lassosum2 <- run_lassosum2(df_beta2)


#### Run C+T ####

run_CT <- function(df_beta) {

  CHR <- as.integer(ukb$map$chromosome)
  POS <- ukb$map$physical.pos
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
      snp = ukb$map$marker.ID,
      A1 = ukb$map$allele1,
      A2 = ukb$map$allele2,
      exclude.ambiguous = FALSE,
      test.bfile = "tmp-data/for_gctb_simu",
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

gctb <- "tmp-data/gctb_2.03beta_Linux/gctb"

corr4gctb <- function(corr) {

  corrT <- as(corr[], "dgTMatrix")

  list_corr <- split(data.frame(i = corrT@i, x = corrT@x),
                     factor(cols_along(corrT))[corrT@j + 1L])

  binfile <- tempfile(tmpdir = "tmp-data", fileext = ".ldm.sparse.bin")
  con <- file(binfile, open = "wb")
  for (df in list_corr) {
    writeBin(as.integer(df$i), con, size = 4)
    writeBin(as.double (df$x), con, size = 4)
  }
  close(con)

  af <- big_colstats(G, ind.row = ind.val, ncores = NCORES)$sum / (2 * length(ind.val))
  info <- ukb$map %>%
    transmute(Chrom = 22, ID = marker.ID, GenPos = 0, PhysPos = physical.pos,
              A1 = allele1, A2 = allele2, A2Freq = af,
              Index = seq_along(list_corr) - 1L,
              WindStart = sapply(list_corr, function(df) df$i[1]),
              WindEnd = sapply(list_corr, function(df) tail(df$i, 1)),
              WindSize = sapply(list_corr, nrow),
              WindWidth = -1,
              N = length(ind.val),
              SamplVar = -1,
              LDsum = Matrix::colSums(corrT))

  bigreadr::fwrite2(info, sub("\\.bin$", ".info", binfile), sep = " ")
  sub("\\.bin$", "", binfile)
}


run_sbayesr <- function(df_beta, corr) {

  # LDSc reg
  corr0 <- readRDS("tmp-data/corr0_simu_val.rds")
  print(ldsc <- snp_ldsc2(corr0, df_beta))
  h2_est <- ldsc[["h2"]]

  af_gwas <- big_colstats(G, ind.row = ind.gwas, ncores = NCORES)$sum / (2 * length(ind.gwas))

  tmp <- tempfile(tmpdir = "tmp-data", fileext = ".ma")
  on.exit(file.remove(tmp), add = TRUE)

  library(dplyr)
  df_beta %>%
    bind_cols(ukb$map) %>%
    transmute(SNP = marker.ID, A1 = allele1, A2 = allele2, freq = af_gwas, b = beta, se = beta_se,
              p = pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE), N = n_eff) %>%
    bigreadr::fwrite2(tmp, sep = " ") %>%
    readLines(n = 5) %>%
    writeLines()

  prefix <- tempfile(tmpdir = "tmp-data")
  on.exit(unlink(paste0(prefix, "*")), add = TRUE)

  if (inherits(corr, "SFBM")) {
    prefix_corr <- corr4gctb(corr)
    on.exit(file.remove(paste0(prefix_corr, c(".bin", ".info"))), add = TRUE)
  } else {
    prefix_corr <- corr
  }

  system(glue::glue(
    gctb,
    " --ldm {prefix_corr}",
    " --gwas-summary {tmp}",
    " --sbayes R",
    " --pi 0.95,0.02,0.02,0.01",
    " --gamma 0.0,0.01,0.1,1",
    " --hsq {h2_est}",
    " --chain-length 5000 --burn-in 1000",
    " --out {prefix} --out-freq 200 --no-mcmc-bin"
  ))  # 40000 SNPs on 1 chromosomes are included.

  tryCatch({
    res_sbayesr <- bigreadr::fread2(paste0(prefix, ".snpRes"))

    stopifnot(all.equal(ukb$map$marker.ID, res_sbayesr$Name))
    stopifnot(all.equal(ukb$map$allele1,   res_sbayesr$A1))
    stopifnot(all.equal(ukb$map$allele2,   res_sbayesr$A2))

    pred_sbayesr <- big_prodVec(G, res_sbayesr$A1Effect, ind.row = ind.test)
    cor(pred_sbayesr, y[ind.test])**2

  }, error = function(e) 0)
}


#### Run PRS-CS-auto ####

run_prscs <- function(df_beta, corr = NULL) {

  tmp <- tempfile(tmpdir = "tmp-data", fileext = ".txt")
  df_beta %>%
    transmute(SNP = ukb$map$marker.ID, A1 = "A", A2 = "C", BETA = beta,
              P = pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE)) %>%
    bigreadr::fwrite2(tmp, sep = "\t") %>%
    readLines(n = 5) %>%
    writeLines()

  on.exit(file.remove(tmp), add = TRUE)

  prefix <- tempfile(tmpdir = "tmp-data")

  dir <- `if`(is.null(corr), "ldblk_ukbb_simu", "ldblk_ukbb_simu_altpop")

  prscs <- "PRScs/PRScs.py"
  system(glue::glue(
    "OMP_NUM_THREADS=", NCORES[[1]],
    " python3 {prscs}",
    " --ref_dir=tmp-data/{dir}",
    " --bim_prefix=tmp-data/{dir}/for_prscs",
    " --sst_file={tmp}",
    " --n_gwas={as.integer(max(df_beta$n_eff))}",
    " --chrom=22",
    " --n_iter=600 --n_burnin=300",
    " --out_dir={prefix}"
  ))

  res_file <- paste0(prefix, "_pst_eff_a1_b0.5_phiauto_chr22.txt")
  on.exit(file.remove(res_file), add = TRUE)

  res <- bigreadr::fread2(res_file)
  stopifnot(identical(res$V2, ukb$map$marker.ID))
  pred_prscs_auto <- big_prodVec(G, res$V6, ind.row = ind.test, ncores = NCORES)
  print(cor(pred_prscs_auto, y[ind.test])**2)
}
