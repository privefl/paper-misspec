library(dplyr)
library(bigsnpr)
ukb <- snp_attach("data/ukbb4simu.rds")
G <- ukb$genotypes


library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results-simu-rounded")
bigassertr::assert_dir("log")


furrr::future_walk(1:10, function(ic) {

  res_file <- paste0("results-simu-rounded/res_", ic, ".rds")
  if (file.exists(res_file)) return(NULL)

  #### Run GWAS and prepare sumstats ####

  load("data/ukbb4simu_ind.RData")

  y <- snp_simuPheno(G, h2 = 0.2, M = 2000, ncores = NCORES)$pheno

  # GWAS to get sumstats
  gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas, ncores = NCORES)
  df_beta <- transmute(gwas, beta = estim, beta_se = std.err, n_eff = length(ind.gwas))

  with(df_beta, min(pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE)))

  # library(ggplot2)
  # maf <- snp_MAF(G, ind.row = ind.gwas, ncores = NCORES)
  # qplot(sqrt(2 * maf * (1 - maf)), 1 / sqrt(n_eff * beta_se^2 + beta^2), alpha = I(0.5),
  #       data = mutate_at(df_beta, c("beta", "beta_se"), ~ signif(., 2))) +
  #   theme_bigstatsr() +
  #   coord_equal() +
  #   scale_color_viridis_d(direction = -1) +
  #   geom_abline(linetype = 2, color = "red") +
  #   labs(x = "Standard deviations derived from allele frequencies",
  #        y = "Standard deviations derived from the summary statistics",
  #        color = "Removed?")


  #### Run them all ####

  source("code/run-methods.R", local = TRUE)

  # Correlation matrices made in 'code/prepare-corr-simu-chr22.R'
  corr1 <- readRDS("tmp-data/corr_simu_val.rds")
  corr2 <- readRDS("tmp-data/corr_simu_val_with_blocks.rds")

  ALL_COR <- list(corr1, corr2)
  ALL_LDPRED <- c("LDpred2-inf", "LDpred2", "LDpred2-low-h2",
                  "LDpred2-auto", "LDpred2-auto-rob")

  all_comb <- tibble::tribble(
    ~ Method,    ~ fun,
    ALL_LDPRED,  run_ldpred2,
    "lassosum2", run_lassosum2,
  ) %>%
    tidyr::expand_grid(corr = ALL_COR, rounded = c(FALSE, TRUE))

  all_comb$r2 <- purrr::pmap(all_comb[-1], function(fun, corr, rounded)  {
    df <- df_beta
    if (rounded) {
      df$beta    <- signif(df$beta,    2)
      df$beta_se <- signif(df$beta_se, 2)
    }
    args <- `if`(is.null(corr), list(df), list(df, corr))
    print(do.call(fun, args))
  })

  saveRDS(select(all_comb, -fun), res_file)
})


library(dplyr)
all_res <- list.files("results-simu-rounded", full.names = TRUE) %>%
  purrr::map_dfr(readRDS) %>%
  mutate(with_blocks = rep(c(FALSE, FALSE, TRUE, TRUE), 18)) %>%
  tidyr::unnest(c(Method, r2)) %>%
  group_by(Method, with_blocks, rounded) %>%
  summarise(r2 = {
    cat(".")
    boot <- replicate(1e3, mean(sample(r2, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(r2), inf = q[[1]], sup = q[[2]]))
  }, .groups = "drop") %>%
  tidyr::unnest_wider("r2") %>%
  print(n = Inf)
#    Method           with_blocks rounded  mean   inf   sup
#  1 lassosum2        FALSE       FALSE   0.173 0.171 0.176
#  2 lassosum2        FALSE       TRUE    0.173 0.170 0.175
#  3 lassosum2        TRUE        FALSE   0.175 0.172 0.177
#  4 lassosum2        TRUE        TRUE    0.174 0.172 0.177
#  5 LDpred2          FALSE       FALSE   0.168 0.166 0.170
#  6 LDpred2          FALSE       TRUE    0.166 0.163 0.169
#  7 LDpred2          TRUE        FALSE   0.175 0.172 0.177
#  8 LDpred2          TRUE        TRUE    0.174 0.171 0.176
#  9 LDpred2-auto     FALSE       FALSE   0.142 0.140 0.144
# 10 LDpred2-auto     FALSE       TRUE    0.142 0.140 0.145
# 11 LDpred2-auto     TRUE        FALSE   0.154 0.151 0.156
# 12 LDpred2-auto     TRUE        TRUE    0.150 0.148 0.152
# 13 LDpred2-auto-rob FALSE       FALSE   0.161 0.159 0.162
# 14 LDpred2-auto-rob FALSE       TRUE    0.160 0.158 0.162
# 15 LDpred2-auto-rob TRUE        FALSE   0.175 0.172 0.177
# 16 LDpred2-auto-rob TRUE        TRUE    0.174 0.172 0.176
# 17 LDpred2-inf      FALSE       FALSE   0.143 0.141 0.146
# 18 LDpred2-inf      FALSE       TRUE    0.143 0.141 0.145
# 19 LDpred2-inf      TRUE        FALSE   0.148 0.145 0.150
# 20 LDpred2-inf      TRUE        TRUE    0.147 0.145 0.150
# 21 LDpred2-low-h2   FALSE       FALSE   0.173 0.171 0.176
# 22 LDpred2-low-h2   FALSE       TRUE    0.173 0.170 0.175
# 23 LDpred2-low-h2   TRUE        FALSE   0.177 0.174 0.179
# 24 LDpred2-low-h2   TRUE        TRUE    0.176 0.173 0.179

ALL_METHODS <- c("lassosum2", "LDpred2", "LDpred2-low-h2",
                 "LDpred2-auto", "LDpred2-auto-rob", "LDpred2-inf")

library(ggplot2)
ggplot(all_res %>% filter(!with_blocks) %>%
         mutate(Method = factor(Method, levels = ALL_METHODS)),
       aes(Method, mean, fill = rounded)) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c("#E69F00", "#56B4E9", "#009E73")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.02), minor_breaks = seq(0, 0.2, by = 0.01),
                     limits = c(0.05, NA), oob = scales::rescale_none) +
  labs(x = "Method", y = "Mean squared correlation between PGS and phenotype",
       fill = "Rounded summary statistics") +
  theme(legend.position = "top") +
  geom_col(data = filter(all_res, with_blocks),
           position = position_dodge(), color = "red",
           alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-rounded.pdf", width = 10, height = 6)
