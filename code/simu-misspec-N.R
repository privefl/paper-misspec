library(dplyr)
library(bigsnpr)
ukb <- snp_attach("data/ukbb4simu.rds")
G <- ukb$genotypes


NCORES <- 14
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results-simu-misN")
bigassertr::assert_dir("log")

res_files <- paste0("results-simu-misN/res_", 1:20, ".rds")

furrr::future_walk(which(!file.exists(res_files)), function(ic) {

  #### Run GWAS and prepare sumstats ####

  load("data/ukbb4simu_ind.RData")

  list_ind <- list(ind.gwas,
                   sort(sample(ind.gwas, length(ind.gwas) * 0.8)),
                   sort(sample(ind.gwas, length(ind.gwas) * 0.6)))

  y <- snp_simuPheno(G, h2 = `if`(ic > 10, 0.04, 0.2), M = 2000, ncores = NCORES)$pheno

  # GWAS to get sumstats
  gwas_set <- sample(rep_len(c(1, 1, 2, 3), ncol(G)))
  df_beta <- data.frame(beta = rep(NA, ncol(G)),
                        beta_se = NA, n_eff = NA, lpval = NA)

  for (k in 1:3) {
    ind_set <- which(gwas_set == k)
    ind.gwas.sub <- list_ind[[k]]
    gwas <- big_univLinReg(G, y[ind.gwas.sub], ind.train = ind.gwas.sub,
                           ind.col = ind_set, ncores = NCORES)
    df_beta$beta[ind_set]    <- gwas$estim
    df_beta$beta_se[ind_set] <- gwas$std.err
    df_beta$n_eff[ind_set]   <- length(ind.gwas.sub)
    df_beta$lpval[ind_set]   <- -predict(gwas)
  }
  with(df_beta, min(pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE)))

  df_beta2 <- df_beta; df_beta2$n_eff <- max(df_beta$n_eff)

  # Imputation of N
  df_beta3 <- df_beta2
  var_G <- big_colstats(G, ind.row = ind.val, ncores = NCORES)$var
  sd_y <- min(with(df_beta2, sqrt(n_eff * beta_se^2 + beta^2))) * sqrt(0.5)
  N_est <- (sd_y^2 / var_G - df_beta2$beta^2) / df_beta2$beta_se^2
  N <- df_beta2$n_eff[1]
  df_beta3$n_eff <- pmin(pmax(0.5 * N, N_est), 1.1 * N)


  #### Run them all ####

  source("code/run-methods.R", local = TRUE)

  ukb <- ukb

  # Correlation matrices made in 'code/prepare-corr-simu-chr22.R'
  corr1 <- readRDS("tmp-data/corr_simu_val.rds")
  corr2 <- readRDS("tmp-data/corr_simu_val_with_blocks.rds")
  corr3 <- readRDS("tmp-data/for_gctb_simu.rds")
  corr4 <- readRDS("tmp-data/for_gctb_simu_with_blocks.rds")

  # run_lassosum(df_beta)
  # run_CT(df_beta)
  # run_prscs(df_beta)
  # run_ldpred2(df_beta, corr3)
  # run_sbayesr(df_beta, corr2)
  # run_sbayesr(df_beta, "tmp-data/for_gctb_simu.ldm.sparse")
  # run_ldpred2(df_beta, corr4)
  # run_sbayesr(df_beta, corr4)

  ALL_N <- c("true", "max", "imputed")
  ALL_LDPRED <- c("LDpred2-inf", "LDpred2", "LDpred2-low-h2",
                  "LDpred2-auto", "LDpred2-auto-rob")

  all_comb <- tibble::tribble(
    ~ Method,      ~ which_N,  ~ fun,         ~ corr,
    "SBayesR",     ALL_N,      run_sbayesr,   list(corr1, corr2, corr3, corr4),
    ALL_LDPRED,    ALL_N,      run_ldpred2,   list(corr1, corr2, corr3, corr4),
    "lassosum2",   ALL_N,      run_lassosum2, list(corr1, corr2, corr3, corr4),
    "lassosum",    ALL_N,      run_lassosum,  list(NULL),
    "C+T",         "any",      run_CT,        list(NULL),
    "PRS-CS-auto", "max",      run_prscs,     list(NULL),
  ) %>%
    tidyr::unnest(which_N) %>%
    tidyr::unnest(corr)

  all_comb$r2 <- purrr::pmap(all_comb[-1], function(which_N, fun, corr)  {
    df <- list(any = df_beta, true = df_beta,
               max = df_beta2, imputed = df_beta3)[[which_N]]
    args <- `if`(is.null(corr), list(df), list(df, corr))
    print(do.call(fun, args))
  })

  saveRDS(select(all_comb, -fun), res_files[ic])
})


ALL_METHODS <- c("C+T", "lassosum", "lassosum2", "LDpred2", "LDpred2-low-h2",
                 "LDpred2-auto", "LDpred2-auto-rob", "LDpred2-inf", "SBayesR", "PRS-CS-auto")

library(dplyr)
all_res <- list.files("results-simu-misN", full.names = TRUE) %>%
  purrr::map_dfr(~ {
    res <- readRDS(.)
    ic <- as.integer(sub(".*res_([0-9]+)\\.rds$", "\\1", .))
    bind_cols(res, h2 = `if`(ic > 10, 0.04, 0.2))
  }) %>%
  rowwise() %>%
  mutate(
    use_blocks = identical(Method, "lassosum") || identical(Method, "PRS-CS-auto") ||
      (!is.null(corr) && grepl("with_blocks", corr$sbk)),
    is_shrunk = !is.null(corr) && grepl("for_gctb", corr$sbk)) %>%
  ungroup() %>%
  tidyr::unnest(c(Method, r2)) %>%
  group_by(Method, which_N, use_blocks, is_shrunk, h2) %>%
  summarise(r2 = {
    cat(".")
    boot <- replicate(1e4, mean(sample(r2, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(r2), inf = q[[1]], sup = q[[2]]))
  }, N = n(), .groups = "drop") %>%
  tidyr::unnest_wider("r2") %>%
  mutate(which_N = factor(which_N, levels = c("any", "true", "max", "imputed")),
         Method  = factor(Method, levels = ALL_METHODS)) %>%
  print(n = Inf)

table(all_res$N)

all_res2 <- filter(all_res, !is_shrunk & h2 == 0.2)

library(ggplot2)
ggplot(filter(all_res2, !use_blocks | Method %in% c("PRS-CS-auto", "lassosum")),
       aes(Method, mean, fill = which_N)) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c(any = "#999999", true = "#E69F00", max = "#56B4E9", imputed = "#009E73")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.02), minor_breaks = seq(0, 0.2, by = 0.01),
                     limits = c(0.05, NA), oob = scales::rescale_none) +
  labs(x = "Method", y = "Mean squared correlation between PGS and phenotype",
       fill = "Sample size") +
  theme(legend.position = "top") +
  geom_col(data = filter(all_res2, use_blocks), position = position_dodge(),
           color = "red", alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-misN.pdf", width = 13, height = 7)

all_res4 <- filter(all_res, !is_shrunk & h2 == 0.04)
ggplot(filter(all_res4, !use_blocks | Method %in% c("PRS-CS-auto", "lassosum")),
       aes(Method, mean, fill = which_N)) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c(any = "#999999", true = "#E69F00", max = "#56B4E9", imputed = "#009E73")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.005), minor_breaks = seq(0, 0.2, by = 0.001)) +
  labs(x = "Method", y = "Mean squared correlation between PGS and phenotype",
       fill = "Sample size") +
  theme(legend.position = "top") +
  geom_col(data = filter(all_res4, use_blocks), position = position_dodge(),
           color = "red", alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-misN-smallh2.pdf", width = 13, height = 7)


all_res3 <- filter(all_res, grepl("SBayesR|LDpred2|lassosum2", Method))
ggplot(filter(all_res3, !use_blocks),
       aes(which_N, mean, fill = ifelse(is_shrunk, "Yes", "No"))) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c("#CC79A7", "#0072B2")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  labs(x = "Sample size", y = "Mean squared correlation between PGS and phenotype",
       fill = "Shrunk LD matrix?") +
  theme(legend.position = "top") +
  facet_grid(h2 ~ Method, scales = "free_y") +
  geom_col(data = filter(all_res3, use_blocks), position = position_dodge(),
           color = "red", alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-misN-shrunk.pdf", width = 13, height = 7)
