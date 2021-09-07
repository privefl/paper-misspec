library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results-simu-info")
bigassertr::assert_dir("log")


furrr::future_walk(1:10, function(ic) {

  res_file <- paste0("results-simu-info/res_", ic, ".rds")
  if (file.exists(res_file)) return(NULL)

  library(dplyr)
  library(bigsnpr)
  simu_true <- snp_attach("data/ukbb4simu.rds")
  G <- simu_true$genotypes

  pheno <- snp_simuPheno(G, h2 = 0.2, M = 2000, ncores = NCORES)
  y <- pheno$pheno

  simu_imp <- snp_attach("data/ukbb4simu_imp.rds")
  INFO <- simu_imp$map$info  # the ones recomputed for the subset
  G_imp <- simu_imp$genotypes

  load("data/ukbb4simu_ind.RData")
  gwas_imp <- big_univLinReg(G_imp, y[ind.gwas], ind.train = ind.gwas,
                             ncores = NCORES)

  df_beta <- transmute(gwas_imp, beta = estim, beta_se = std.err,
                       n_eff = length(ind.gwas))

  df_beta2 <- mutate(df_beta,
                     beta = beta * sqrt(INFO),
                     beta_se = beta_se * sqrt(INFO))
  df_beta3 <- mutate(df_beta,
                     beta = beta * INFO,
                     n_eff = n_eff * INFO)
  df_beta4 <- mutate(df_beta,
                     beta = beta * INFO,
                     beta_se = beta_se * sqrt(INFO),
                     n_eff = n_eff * INFO)


  #### Run them all ####

  source("code/run-methods.R", local = TRUE)

  # Correlation matrices made in 'code/prepare-corr-simu-chr22.R'
  corr0 <- readRDS("tmp-data/corr0_simu_val.rds")
  corr1 <- readRDS("tmp-data/corr_simu_val.rds")
  corr2 <- readRDS("tmp-data/corr_simu_val_with_blocks.rds")
  dim(corr0)  # 40000 x 40000

  ALL_LDPRED <- c("LDpred2-inf", "LDpred2", "LDpred2-low-h2",
                  "LDpred2-auto", "LDpred2-auto-rob")

  all_comb <- tibble::tribble(
    ~ Method,     ~ fun,
    ALL_LDPRED,   run_ldpred2,
    "lassosum2",  run_lassosum2,
  ) %>%
    tidyr::expand_grid(correction = c("none", "sqrt_info", "info", "in_between"),
                       corr = list(corr1, corr2))

  all_comb$r2 <- purrr::pmap(all_comb[-1], function(correction, fun, corr)  {
    df <- list(none = df_beta, sqrt_info = df_beta2,
               info = df_beta3, in_between = df_beta4)[[correction]]
    print(fun(df, corr))
  })

  saveRDS(select(all_comb, -fun), res_file)
})


library(dplyr)
all_res <- list.files("results-simu-info", full.names = TRUE) %>%
  purrr::map_dfr(readRDS) %>%
  rowwise() %>%
  mutate(corr = ifelse(is.null(corr), "none",
                       ifelse(grepl("with_blocks", corr$sbk), "with blocks", "normal"))) %>%
  ungroup() %>%
  tidyr::unnest(c(Method, r2)) %>%
  group_by(Method, correction, corr) %>%
  summarise(r2 = {
    cat(".")
    boot <- replicate(1e4, mean(sample(r2, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(r2), inf = q[[1]], sup = q[[2]]))
  }, .groups = "drop") %>%
  tidyr::unnest_wider("r2") %>%
  mutate(correction = factor(correction,
                             levels = c("none", "sqrt_info", "info", "in_between"))) %>%
  print(n = Inf)
#    Method           correction corr         mean   inf   sup
#  1 lassosum2        in_between normal      0.157 0.155 0.159
#  2 lassosum2        in_between with blocks 0.158 0.157 0.160
#  3 lassosum2        info       normal      0.157 0.156 0.159
#  4 lassosum2        info       with blocks 0.159 0.157 0.160
#  5 lassosum2        none       normal      0.154 0.152 0.156
#  6 lassosum2        none       with blocks 0.156 0.154 0.157
#  7 lassosum2        sqrt_info  normal      0.157 0.155 0.159
#  8 lassosum2        sqrt_info  with blocks 0.158 0.157 0.160
#  9 LDpred2          in_between normal      0.155 0.153 0.157
# 10 LDpred2          in_between with blocks 0.159 0.158 0.160
# 11 LDpred2          info       normal      0.155 0.153 0.157
# 12 LDpred2          info       with blocks 0.158 0.156 0.159
# 13 LDpred2          none       normal      0.150 0.149 0.152
# 14 LDpred2          none       with blocks 0.156 0.155 0.158
# 15 LDpred2          sqrt_info  normal      0.153 0.151 0.154
# 16 LDpred2          sqrt_info  with blocks 0.159 0.158 0.161
# 17 LDpred2-auto     in_between normal      0.136 0.134 0.137
# 18 LDpred2-auto     in_between with blocks 0.158 0.157 0.160
# 19 LDpred2-auto     info       normal      0.154 0.152 0.155
# 20 LDpred2-auto     info       with blocks 0.159 0.157 0.160
# 21 LDpred2-auto     none       normal      0.129 0.128 0.131
# 22 LDpred2-auto     none       with blocks 0.150 0.148 0.152
# 23 LDpred2-auto     sqrt_info  normal      0.133 0.131 0.135
# 24 LDpred2-auto     sqrt_info  with blocks 0.153 0.151 0.155
# 25 LDpred2-auto-rob in_between normal      0.156 0.155 0.158
# 26 LDpred2-auto-rob in_between with blocks 0.160 0.159 0.162
# 27 LDpred2-auto-rob info       normal      0.157 0.156 0.159
# 28 LDpred2-auto-rob info       with blocks 0.159 0.157 0.160
# 29 LDpred2-auto-rob none       normal      0.145 0.143 0.147
# 30 LDpred2-auto-rob none       with blocks 0.155 0.154 0.157
# 31 LDpred2-auto-rob sqrt_info  normal      0.149 0.147 0.151
# 32 LDpred2-auto-rob sqrt_info  with blocks 0.159 0.157 0.160
# 33 LDpred2-inf      in_between normal      0.136 0.134 0.137
# 34 LDpred2-inf      in_between with blocks 0.138 0.137 0.140
# 35 LDpred2-inf      info       normal      0.137 0.135 0.139
# 36 LDpred2-inf      info       with blocks 0.139 0.138 0.141
# 37 LDpred2-inf      none       normal      0.130 0.128 0.131
# 38 LDpred2-inf      none       with blocks 0.134 0.133 0.136
# 39 LDpred2-inf      sqrt_info  normal      0.133 0.132 0.135
# 40 LDpred2-inf      sqrt_info  with blocks 0.137 0.136 0.139
# 41 LDpred2-low-h2   in_between normal      0.158 0.156 0.160
# 42 LDpred2-low-h2   in_between with blocks 0.160 0.159 0.161
# 43 LDpred2-low-h2   info       normal      0.156 0.155 0.158
# 44 LDpred2-low-h2   info       with blocks 0.158 0.156 0.159
# 45 LDpred2-low-h2   none       normal      0.155 0.153 0.157
# 46 LDpred2-low-h2   none       with blocks 0.158 0.156 0.160
# 47 LDpred2-low-h2   sqrt_info  normal      0.157 0.156 0.159
# 48 LDpred2-low-h2   sqrt_info  with blocks 0.160 0.159 0.162

ALL_METHODS <- c("C+T", "lassosum", "lassosum2", "LDpred2", "LDpred2-low-h2",
                 "LDpred2-auto", "LDpred2-auto-rob", "LDpred2-inf")

library(ggplot2)
ggplot(all_res %>% filter(corr != "with blocks") %>%
         mutate(Method = factor(Method, levels = ALL_METHODS)),
       aes(Method, mean, fill = as.factor(correction))) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c("#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.02), minor_breaks = seq(0, 0.2, by = 0.01),
                     limits = c(0.08, NA), oob = scales::rescale_none) +
  labs(x = "Method", y = "Mean squared correlation between PGS and phenotype",
       fill = "Correction") +
  theme(legend.position = "top") +
  geom_col(data = filter(all_res, corr == "with blocks"),
           position = position_dodge(), color = "red",
           alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-info.pdf", width = 10, height = 6)
