library(dplyr)
library(bigsnpr)
ukb <- snp_attach("data/ukbb4simu.rds")
G <- ukb$genotypes


library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results-simu-altpop")
bigassertr::assert_dir("log")


furrr::future_walk(1:10, function(ic) {

  res_file <- paste0("results-simu-altpop/res_", ic, ".rds")
  if (file.exists(res_file)) return(NULL)

  #### Run GWAS and prepare sumstats ####

  load("data/ukbb4simu_ind.RData")

  y <- snp_simuPheno(G, h2 = 0.2, M = 2000, ncores = NCORES)$pheno

  # GWAS to get sumstats
  gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas, ncores = NCORES)
  df_beta <- transmute(gwas, beta = estim, beta_se = std.err, n_eff = length(ind.gwas))

  with(df_beta, min(pchisq((beta / beta_se)^2, df = 1, lower.tail = FALSE)))

  #### Run them all ####

  source("code/run-methods.R", local = TRUE)

  # Correlation matrices made in 'code/prepare-corr-simu-chr22.R'
  corr1 <- readRDS("tmp-data/corr_simu_val.rds")
  corr2 <- readRDS("tmp-data/corr_simu_val_with_blocks.rds")
  # Correlation matrices made in 'code/prepare-corr-altpop.R'
  corr3 <- readRDS("tmp-data/corr_simu_altpop.rds")
  corr4 <- readRDS("tmp-data/corr_simu_altpop_with_blocks.rds")

  ALL_COR <- list(corr1, corr2, corr3, corr4)
  ALL_LDPRED <- c("LDpred2-inf", "LDpred2", "LDpred2-low-h2",
                  "LDpred2-auto", "LDpred2-auto-rob")

  all_comb <- tibble::tribble(
    ~ Method,    ~ fun,         ~ corr,
    ALL_LDPRED,  run_ldpred2,   ALL_COR,
    "lassosum2", run_lassosum2, ALL_COR
  ) %>%
    tidyr::unnest(corr)

  all_comb$r2 <- purrr::pmap(all_comb[-1], function(fun, corr)  {
    df <- df_beta
    args <- `if`(is.null(corr), list(df), list(df, corr))
    print(do.call(fun, args))
  })

  saveRDS(select(all_comb, -fun), res_file)
})


library(dplyr)
all_res <- list.files("results-simu-altpop", full.names = TRUE) %>%
  purrr::map_dfr(readRDS) %>%
  mutate(with_blocks = rep(c(FALSE, TRUE, FALSE, TRUE), 20),
         altpop = rep(c("No", "No", "Yes", "Yes"), 20)) %>%
  tidyr::unnest(c(Method, r2)) %>%
  group_by(Method, with_blocks, altpop) %>%
  summarise(r2 = {
    cat(".")
    boot <- replicate(1e4, mean(sample(r2, replace = TRUE)))
    q <- quantile(boot, c(0.025, 0.975))
    list(c(mean = mean(r2), inf = q[[1]], sup = q[[2]]))
  }, .groups = "drop") %>%
  tidyr::unnest_wider("r2") %>%
  print(n = Inf)
#    Method           with_blocks altpop  mean   inf   sup
#  1 lassosum2        FALSE       No     0.174 0.171 0.177
#  2 lassosum2        FALSE       Yes    0.167 0.164 0.170
#  3 lassosum2        TRUE        No     0.175 0.173 0.178
#  4 lassosum2        TRUE        Yes    0.168 0.165 0.171
#  5 LDpred2          FALSE       No     0.168 0.165 0.171
#  6 LDpred2          FALSE       Yes    0.147 0.145 0.149
#  7 LDpred2          TRUE        No     0.175 0.172 0.177
#  8 LDpred2          TRUE        Yes    0.151 0.149 0.154
#  9 LDpred2-auto     FALSE       No     0.143 0.141 0.146
# 10 LDpred2-auto     FALSE       Yes    0.139 0.137 0.141
# 11 LDpred2-auto     TRUE        No     0.153 0.150 0.155
# 12 LDpred2-auto     TRUE        Yes    0.142 0.140 0.144
# 13 LDpred2-auto-rob FALSE       No     0.162 0.159 0.166
# 14 LDpred2-auto-rob FALSE       Yes    0.146 0.143 0.148
# 15 LDpred2-auto-rob TRUE        No     0.175 0.172 0.178
# 16 LDpred2-auto-rob TRUE        Yes    0.151 0.149 0.154
# 17 LDpred2-inf      FALSE       No     0.144 0.142 0.146
# 18 LDpred2-inf      FALSE       Yes    0.140 0.138 0.142
# 19 LDpred2-inf      TRUE        No     0.148 0.146 0.151
# 20 LDpred2-inf      TRUE        Yes    0.142 0.140 0.144
# 21 LDpred2-low-h2   FALSE       No     0.174 0.171 0.177
# 22 LDpred2-low-h2   FALSE       Yes    0.169 0.167 0.172
# 23 LDpred2-low-h2   TRUE        No     0.177 0.175 0.180
# 24 LDpred2-low-h2   TRUE        Yes    0.170 0.168 0.173

ALL_METHODS <- c("lassosum2", "LDpred2", "LDpred2-low-h2",
                 "LDpred2-auto", "LDpred2-auto-rob", "LDpred2-inf")

library(ggplot2)
ggplot(all_res %>% filter(!with_blocks) %>%
         mutate(Method = factor(Method, levels = ALL_METHODS)),
       aes(Method, mean, fill = altpop)) +
  bigstatsr::theme_bigstatsr(0.8) +
  scale_fill_manual(values =  c("#E69F00", "#56B4E9")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = inf, ymax = sup),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.02), minor_breaks = seq(0, 0.2, by = 0.01),
                     limits = c(0.05, NA), oob = scales::rescale_none) +
  labs(x = "Method", y = "Mean squared correlation between PGS and phenotype",
       fill = "Alternative LD reference") +
  theme(legend.position = "top") +
  geom_col(data = filter(all_res, with_blocks),
           position = position_dodge(), color = "red",
           alpha = 0, show.legend = FALSE)
# ggsave("figures/simu-altpop.pdf", width = 10, height = 6)
