library(bigsnpr)
ukb <- snp_attach("data/ukbb4simu.rds")
G <- ukb$genotypes

load("data/ukbb4simu_ind.RData")

(NCORES <- nb_cores())

# Correlation matrix made in 'code/prepare-corr-simu-chr22.R'
corr <- readRDS("tmp-data/corr_simu_val.rds")

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
params$score <- apply(pred_lassosum[ind.test, ], 2, cor, y = y[ind.test])

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

plot_grid(
  qplot(auto_score, score, color = as.factor(delta), data = params) +
    theme_bw(15) +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Score from pseudo-validation (using correlations)",
         y = "Score from validation", color = expression(delta)) +
    theme(legend.position = c(0.15, 0.7)),
  qplot(auto_score2, score, color = as.factor(delta), data = params) +
    theme_bw(15) +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Score from pseudo-validation (using p-values)",
         y = "Score from validation", color = expression(delta)) +
    theme(legend.position = c(0.8, 0.3)),
  scale = 0.95, labels = c("A", "B"), label_size = 16, ncol = 1
)
# ggsave("figures/pseudoval.pdf", width = 7, height = 10)
