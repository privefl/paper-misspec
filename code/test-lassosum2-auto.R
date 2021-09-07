library(bigsnpr)
ukb <- snp_attach("data/ukbb4simu.rds")
G <- ukb$genotypes

load("data/ukbb4simu_ind.RData")

(NCORES <- nb_cores())

corr0 <- runonce::save_run({
  POS2 <- snp_asGeneticPos(as.integer(ukb$map$chromosome),
                           ukb$map$physical.pos, dir = "tmp-data")
  snp_cor(G, ind.row = ind.val, infos.pos = POS2, size = 3 / 1000, ncores = nb_cores())
}, file = "tmp-data/corr0_simu_val.rds")
corr <- runonce::save_run(as_SFBM(corr0, "tmp-data/corr_simu_val"),
                          file = "tmp-data/corr_simu_val.rds")

y <- snp_simuPheno(G, h2 = 0.2, M = 2000, ncores = NCORES)$pheno

# GWAS to get sumstats
gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas, ncores = NCORES)
library(dplyr)
df_beta <- gwas %>%
  transmute(beta = estim, beta_se = std.err, n_eff = length(ind.gwas))

# lassosum2
beta_lassosum <- snp_lassosum2(corr, df_beta, ncores = NCORES)
(params <- attr(beta_lassosum, "grid_param"))

# validation
pred_lassosum <- big_prodMat(G, beta_lassosum, ncores = NCORES)
params$score <- apply(pred_lassosum[ind.val, ], 2, cor, y = y[ind.val])

# pseudo-validation
scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
beta_hat <- df_beta$beta / scale

fdr <- fdrtool::fdrtool(beta_hat, statistic = "correlation", plot = FALSE)
beta_hat_shrunk <- round(beta_hat * (1 - fdr$lfdr), 16)

params$auto_score <- apply(beta_lassosum, 2, function(beta) {
  cat(".")
  beta <- beta / scale
  bRb <- crossprod(beta, bigsparser::sp_prodVec(corr, beta))
  crossprod(beta, beta_hat_shrunk) / sqrt(bRb)
})

library(ggplot2)
qplot(auto_score, score, color = delta, data = params) +
  theme_bw(15) +
  scale_color_viridis_c(trans = "log10") +
  labs(x = "Score from pseudo-validation", y = "Score from validation")


pval <- predict(gwas, log10 = FALSE)
fdr2 <- fdrtool::fdrtool(pval, statistic = "pvalue", plot = FALSE)
beta_hat_shrunk2 <- beta_hat * (1 - fdr2$lfdr)

params$auto_score2 <- apply(beta_lassosum, 2, function(beta) {
  cat(".")
  beta <- beta / scale
  bRb <- crossprod(beta, bigsparser::sp_prodVec(corr, beta))
  crossprod(beta, beta_hat_shrunk2) / sqrt(bRb)
})

source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
plot_grid2(list(
  qplot(auto_score, score, color = as.factor(delta), data = params) +
    theme_bw(15) +
    labs(x = "Score from pseudo-validation (using correlations)",
         y = "Score from validation", color = "delta"),
  qplot(auto_score2, score, color = as.factor(delta), data = params) +
    theme_bw(15) +
    labs(x = "Score from pseudo-validation (using p-values)",
         y = "Score from validation")
), scale = 0.95, labels = c("A", "B"), label_size = 16, ncol = 1,
title_ratio = 0, legend_ratio = 0.15)
# ggsave("figures/pseudoval.pdf", width = 8, height = 10)
