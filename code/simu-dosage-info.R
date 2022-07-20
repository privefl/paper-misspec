library(future.batchtools)
NCORES <- 14
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results-simu-info")
bigassertr::assert_dir("log")

res_files <- paste0("results-simu-info/res_", 1:20, ".rds")

furrr::future_walk(which(!file.exists(res_files)), function(ic) {

  library(dplyr)
  library(bigsnpr)
  ukb <- simu_true <- snp_attach("data/ukbb4simu.rds")
  G <- simu_true$genotypes

  pheno <- snp_simuPheno(G, h2 = `if`(ic > 10, 0.04, 0.2), M = 2000, ncores = NCORES)
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
  corr1 <- readRDS("tmp-data/corr_simu_val.rds")
  corr2 <- readRDS("tmp-data/corr_simu_val_with_blocks.rds")

  ALL_LDPRED <- c("LDpred2-inf", "LDpred2", "LDpred2-low-h2",
                  "LDpred2-auto", "LDpred2-auto-rob")

  all_comb <- tibble::tribble(
    ~ Method,      ~ fun,         ~ corr,
    "PRS-CS-auto", run_prscs,     list(NULL),
    ALL_LDPRED,    run_ldpred2,   list(corr1, corr2),
    "lassosum2",   run_lassosum2, list(corr1, corr2),
  ) %>%
    tidyr::unnest(corr) %>%
    tidyr::expand_grid(correction = c("none", "sqrt_info", "info", "in_between"))

  all_comb$r2 <- purrr::pmap(all_comb[-1], function(correction, fun, corr)  {
    df <- list(none = df_beta, sqrt_info = df_beta2,
               info = df_beta3, in_between = df_beta4)[[correction]]
    print(fun(df, corr))
  })

  saveRDS(select(all_comb, -fun), res_files[ic])
})


ALL_METHODS <- c("lassosum2", "LDpred2", "LDpred2-low-h2",
                 "LDpred2-auto", "LDpred2-auto-rob", "LDpred2-inf", "PRS-CS-auto")

library(dplyr)
all_res <- list.files("results-simu-info", full.names = TRUE) %>%
  purrr::map_dfr(~ {
    res <- readRDS(.)
    ic <- as.integer(sub(".*res_([0-9]+)\\.rds$", "\\1", .))
    bind_cols(res, h2 = `if`(ic > 10, 0.04, 0.2))
  }) %>%
  rowwise() %>%
  mutate(use_blocks = identical(Method, "PRS-CS-auto") || grepl("with_blocks", corr$sbk)) %>%
  ungroup() %>%
  tidyr::unnest(c(Method, r2)) %>%
  group_by(Method, correction, use_blocks, h2) %>%
  summarise(r2 = {
    cat(".")
    boot <- replicate(1e4, mean(sample(r2, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(r2), inf = q[[1]], sup = q[[2]]))
  }, N = n(), .groups = "drop") %>%
  tidyr::unnest_wider("r2") %>%
  mutate(correction = factor(correction, levels = c("none", "sqrt_info", "info", "in_between")),
         Method  = factor(Method, levels = ALL_METHODS)) %>%
  print(n = Inf)

table(all_res$N)

all_res2 <- filter(all_res, h2 == 0.2)

library(ggplot2)
ggplot(filter(all_res2, !use_blocks | Method == "PRS-CS-auto"),
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
  geom_col(data = filter(all_res2, use_blocks), position = position_dodge(),
           color = "red", alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-info.pdf", width = 10, height = 6)

all_res3 <- filter(all_res, h2 == 0.04)
ggplot(filter(all_res3, !use_blocks | Method == "PRS-CS-auto"),
       aes(Method, mean, fill = as.factor(correction))) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c("#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.005), minor_breaks = seq(0, 0.2, by = 0.001)) +
  labs(x = "Method", y = "Mean squared correlation between PGS and phenotype",
       fill = "Correction") +
  theme(legend.position = "top") +
  geom_col(data = filter(all_res3, use_blocks), position = position_dodge(),
           color = "red", alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-info-smallh2.pdf", width = 10, height = 6)
