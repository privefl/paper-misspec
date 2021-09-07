library(dplyr)
library(bigsnpr)
ukb <- snp_attach("data/ukbb4simu.rds")
G <- ukb$genotypes

load("data/ukbb4simu_ind.RData")

# prepare data in bed format
data <- snp_fake(n = nrow(G), m = 1)
data$map <- ukb$map
data$map$genetic.dist <- snp_asGeneticPos(as.integer(ukb$map$chromosome),
                                          ukb$map$physical.pos, dir = "tmp-data")
data$map$marker.ID <- paste0("SNP", rows_along(data$map))
data$genotypes <- G
# snp_writeBed(data, ind.row = ind.val, bedfile = "tmp-data/simu_chr22.bed")

library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results-simu-misN")
bigassertr::assert_dir("log")


furrr::future_walk(1:10, function(ic) {

  res_file <- paste0("results-simu-misN/res_", ic, ".rds")
  if (file.exists(res_file)) return(NULL)

  #### Run GWAS and prepare sumstats ####

  load("data/ukbb4simu_ind.RData")

  list_ind <- list(ind.gwas,
                   sort(sample(ind.gwas, length(ind.gwas) * 0.8)),
                   sort(sample(ind.gwas, length(ind.gwas) * 0.6)))

  y <- snp_simuPheno(G, h2 = 0.2, M = 2000, ncores = NCORES)$pheno

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

  data <- data

  # Correlation matrices made in 'code/prepare-corr-simu-chr22.R'
  corr0 <- readRDS("tmp-data/corr0_simu_val.rds")
  corr1 <- readRDS("tmp-data/corr_simu_val.rds")
  corr2 <- readRDS("tmp-data/corr_simu_val_with_blocks.rds")
  dim(corr0)  # 40000 x 40000

  ALL_N <- c("true", "max", "imputed")
  ALL_LDPRED <- c("LDpred2-inf", "LDpred2", "LDpred2-low-h2",
                  "LDpred2-auto", "LDpred2-auto-rob")

  all_comb <- tibble::tribble(
    ~ Method,    ~ which_N,  ~ fun,         ~ corr,
    ALL_LDPRED,  ALL_N,      run_ldpred2,   list(corr1, corr2),
    "lassosum2", ALL_N,      run_lassosum2, list(corr1, corr2),
    "lassosum",  ALL_N,      run_lassosum,  list(NULL),
    "C+T",       "any",      run_CT,        list(NULL),
  ) %>%
    tidyr::unnest(which_N) %>%
    tidyr::unnest(corr)

  all_comb$r2 <- purrr::pmap(all_comb[-1], function(which_N, fun, corr)  {
    df <- list(any = df_beta, true = df_beta,
               max = df_beta2, imputed = df_beta3)[[which_N]]
    args <- `if`(is.null(corr), list(df), list(df, corr))
    print(do.call(fun, args))
  })

  saveRDS(select(all_comb, -fun), res_file)
})


library(dplyr)
all_res <- list.files("results-simu-misN", full.names = TRUE) %>%
  purrr::map_dfr(readRDS) %>%
  rowwise() %>%
  mutate(corr = ifelse(is.null(corr), "none",
                       ifelse(grepl("with_blocks", corr$sbk), "with blocks", "normal"))) %>%
  ungroup() %>%
  tidyr::unnest(c(Method, r2)) %>%
  group_by(Method, which_N, corr) %>%
  summarise(r2 = {
    cat(".")
    boot <- replicate(1e4, mean(sample(r2, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(r2), inf = q[[1]], sup = q[[2]]))
  }, .groups = "drop") %>%
  tidyr::unnest_wider("r2") %>%
  mutate(which_N = factor(which_N, levels = c("any", "true", "max", "imputed"))) %>%
  print(n = Inf)
#    Method           which_N corr         mean   inf   sup
#  1 C+T              any     none        0.123 0.121  0.125
#  2 lassosum         imputed none        0.161 0.159  0.163
#  3 lassosum         max     none        0.157 0.155  0.160
#  4 lassosum         true    none        0.161 0.159  0.163
#  5 lassosum2        imputed normal      0.169 0.167  0.170
#  6 lassosum2        imputed with blocks 0.170 0.168  0.172
#  7 lassosum2        max     normal      0.163 0.159  0.166
#  8 lassosum2        max     with blocks 0.164 0.160  0.167
#  9 lassosum2        true    normal      0.169 0.167  0.170
# 10 lassosum2        true    with blocks 0.170 0.168  0.172
# 11 LDpred2          imputed normal      0.159 0.157  0.161
# 12 LDpred2          imputed with blocks 0.163 0.161  0.165
# 13 LDpred2          max     normal      0.134 0.130  0.137
# 14 LDpred2          max     with blocks 0.136 0.132  0.140
# 15 LDpred2          true    normal      0.159 0.157  0.161
# 16 LDpred2          true    with blocks 0.163 0.161  0.165
# 17 LDpred2-auto     imputed normal      0.140 0.138  0.142
# 18 LDpred2-auto     imputed with blocks 0.143 0.141  0.145
# 19 LDpred2-auto     max     normal      0.119 0.111  0.124
# 20 LDpred2-auto     max     with blocks 0.120 0.111  0.126
# 21 LDpred2-auto     true    normal      0.140 0.138  0.142
# 22 LDpred2-auto     true    with blocks 0.143 0.142  0.145
# 23 LDpred2-auto-rob imputed normal      0.151 0.149  0.154
# 24 LDpred2-auto-rob imputed with blocks 0.165 0.163  0.167
# 25 LDpred2-auto-rob max     normal      0.109 0.0932 0.120
# 26 LDpred2-auto-rob max     with blocks 0.111 0.0925 0.125
# 27 LDpred2-auto-rob true    normal      0.152 0.150  0.154
# 28 LDpred2-auto-rob true    with blocks 0.165 0.163  0.167
# 29 LDpred2-inf      imputed normal      0.141 0.139  0.143
# 30 LDpred2-inf      imputed with blocks 0.144 0.142  0.145
# 31 LDpred2-inf      max     normal      0.123 0.117  0.127
# 32 LDpred2-inf      max     with blocks 0.126 0.120  0.130
# 33 LDpred2-inf      true    normal      0.141 0.139  0.143
# 34 LDpred2-inf      true    with blocks 0.144 0.142  0.146
# 35 LDpred2-low-h2   imputed normal      0.169 0.167  0.171
# 36 LDpred2-low-h2   imputed with blocks 0.172 0.170  0.174
# 37 LDpred2-low-h2   max     normal      0.163 0.159  0.166
# 38 LDpred2-low-h2   max     with blocks 0.163 0.160  0.166
# 39 LDpred2-low-h2   true    normal      0.169 0.167  0.171
# 40 LDpred2-low-h2   true    with blocks 0.172 0.171  0.174
ALL_METHODS <- c("C+T", "lassosum", "lassosum2", "LDpred2", "LDpred2-low-h2",
                 "LDpred2-auto", "LDpred2-auto-rob", "LDpred2-inf")

library(ggplot2)
ggplot(all_res %>% filter(corr != "with blocks") %>%
         mutate(Method = factor(Method, levels = ALL_METHODS)),
       aes(Method, mean, fill = which_N)) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c(any = "#999999", "#E69F00", "#56B4E9", "#009E73")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.02), minor_breaks = seq(0, 0.2, by = 0.01),
                     limits = c(0.05, NA), oob = scales::rescale_none) +
  labs(x = "Method", y = "Mean squared correlation between PGS and phenotype",
       fill = "Sample size") +
  theme(legend.position = "top") +
  geom_col(data = filter(all_res, corr == "with blocks"),
           position = position_dodge(), color = "red",
           alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-misN.pdf", width = 11, height = 6)
