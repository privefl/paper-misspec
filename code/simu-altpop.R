library(dplyr)
library(bigsnpr)
ukb <- snp_attach("data/ukbb4simu.rds")
G <- ukb$genotypes


library(future.batchtools)
NCORES <- 14
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results-simu-altpop")
bigassertr::assert_dir("log")

res_files <- paste0("results-simu-altpop/res_", 1:20, ".rds")

furrr::future_walk(which(!file.exists(res_files)), function(ic) {

  #### Run GWAS and prepare sumstats ####

  load("data/ukbb4simu_ind.RData")

  y <- snp_simuPheno(G, h2 = `if`(ic > 10, 0.04, 0.2), M = 2000, ncores = NCORES)$pheno

  # GWAS to get sumstats
  gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas, ncores = NCORES)
  df_beta <- transmute(gwas, beta = estim, beta_se = std.err, n_eff = length(ind.gwas))

  with(df_beta, min(pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE)))

  #### Run them all ####

  source("code/run-methods.R", local = TRUE)

  ukb <- ukb

  # Correlation matrices made in 'code/prepare-corr-simu-chr22.R'
  corr1 <- readRDS("tmp-data/corr_simu_val.rds")
  corr2 <- readRDS("tmp-data/corr_simu_val_with_blocks.rds")
  # Correlation matrices made in 'code/prepare-corr-altpop.R'
  corr3 <- readRDS("tmp-data/corr_simu_altpop.rds")
  corr4 <- readRDS("tmp-data/corr_simu_altpop_with_blocks.rds")
  corr5 <- readRDS("tmp-data/for_gctb_simu_altpop.rds")
  corr6 <- readRDS("tmp-data/for_gctb_simu_altpop_with_blocks.rds")

  ALL_COR <- list(corr1, corr2, corr3, corr4, corr5, corr6)
  ALL_LDPRED <- c("LDpred2-inf", "LDpred2", "LDpred2-low-h2",
                  "LDpred2-auto", "LDpred2-auto-rob")

  all_comb <- tibble::tribble(
    ~ Method,      ~ fun,         ~ corr,
    ALL_LDPRED,    run_ldpred2,   ALL_COR,
    "lassosum2",   run_lassosum2, ALL_COR,
    "PRS-CS-auto", run_prscs,     list(NULL, "altpop"),
  ) %>%
    tidyr::unnest(corr)

  all_comb$r2 <- purrr::pmap(all_comb[-1], function(fun, corr)  {
    df <- df_beta
    args <- `if`(is.null(corr), list(df), list(df, corr))
    print(do.call(fun, args))
  })

  saveRDS(select(all_comb, -fun), res_files[ic])
})


ALL_METHODS <- c("lassosum2", "LDpred2", "LDpred2-low-h2",
                 "LDpred2-auto", "LDpred2-auto-rob", "LDpred2-inf", "PRS-CS-auto")

library(dplyr)
all_res <- list.files("results-simu-altpop", full.names = TRUE) %>%
  purrr::map_dfr(~ {
    res <- readRDS(.)
    ic <- as.integer(sub(".*res_([0-9]+)\\.rds$", "\\1", .))
    bind_cols(res, h2 = `if`(ic > 10, 0.04, 0.2))
  }) %>%
  rowwise() %>%
  mutate(
    use_blocks = Method == "PRS-CS-auto" || grepl("with_blocks", corr$sbk),
    corr = {
      if (inherits(corr, "SFBM")) {
        `if`(grepl("for_gctb", corr$sbk), "S. Eur (shrunk)",
             `if`(grepl("altpop", corr$sbk), "S. Eur", "N.W. Eur"))
      } else {
        `if`(is.null(corr), "N.W. Eur", "S. Eur")
      }
    }) %>%
  ungroup() %>%
  tidyr::unnest(c(Method, r2)) %>%
  mutate(Method = factor(Method, levels = ALL_METHODS)) %>%
  group_by(Method, use_blocks, corr, h2) %>%
  summarise(r2 = {
    cat(".")
    boot <- replicate(1e4, mean(sample(r2, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(r2), inf = q[[1]], sup = q[[2]]))
  }, N = n(), .groups = "drop") %>%
  tidyr::unnest_wider("r2") %>%
  print(n = Inf)

table(all_res$N)

all_res2 <- filter(all_res, h2 == 0.2)

library(ggplot2)
ggplot(filter(all_res2, !use_blocks | Method == "PRS-CS-auto"), aes(Method, mean, fill = corr)) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c("#999999", "#009E73", "#F0E442")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.02), minor_breaks = seq(0, 0.2, by = 0.01),
                     limits = c(0.05, NA), oob = scales::rescale_none) +
  labs(x = "Method", y = "Mean squared correlation between PGS and phenotype",
       fill = "LD reference") +
  theme(legend.position = "top") +
  geom_col(data = filter(all_res2, use_blocks), position = position_dodge(),
           color = "red", alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-altpop.pdf", width = 10, height = 6)


all_res3 <- filter(all_res, h2 == 0.04)
ggplot(filter(all_res3, !use_blocks | Method == "PRS-CS-auto"), aes(Method, mean, fill = corr)) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c("#999999", "#009E73", "#F0E442")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.005), minor_breaks = seq(0, 0.2, by = 0.001)) +
  labs(x = "Method", y = "Mean squared correlation between PGS and phenotype",
       fill = "LD reference") +
  theme(legend.position = "top") +
  geom_col(data = filter(all_res3, use_blocks), position = position_dodge(),
           color = "red", alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-altpop-smallh2.pdf", width = 10, height = 6)
