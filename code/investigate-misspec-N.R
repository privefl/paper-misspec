library(dplyr)
library(ggplot2)
library(bigsnpr)
G <- snp_attach("data/ukbb4simu.rds")$genotypes

load("data/ukbb4simu_ind.RData")

# Correlation matrices made in 'code/prepare-corr-simu-chr22.R'
corr0 <- readRDS("tmp-data/corr0_simu_val.rds")
corr1 <- readRDS("tmp-data/corr_simu_val.rds")
corr2 <- readRDS("tmp-data/corr_simu_val_with_blocks.rds")
dim(corr0)  # 40000 x 40000

NCORES <- nb_cores()


#### Run GWAS and prepare sumstats ####

set.seed(1)
y <- snp_simuPheno(G, h2 = 0.2, M = 2000)$pheno

list_ind <- list(ind.gwas,
                 sort(sample(ind.gwas, length(ind.gwas) * 0.8)),
                 sort(sample(ind.gwas, length(ind.gwas) * 0.6)))

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

df_beta2 <- df_beta; df_beta2$n_eff <- max(df_beta$n_eff)  # using max(N)

# Quality control plot from the LDpred2 paper
sd_val <- sqrt(big_colstats(G, ind.row = ind.val, ncores = NCORES)$var)
sd_ss <- 1 / with(df_beta2, sqrt(n_eff * beta_se^2 + beta^2))  # using max(N)
sd_y <- sqrt(0.5) / max(sd_ss)
sd_ss <- sd_y * sd_ss

qplot(sd_val, sd_ss, alpha = I(0.6),
      color = factor(df_beta$n_eff, levels = sort(unique(df_beta$n_eff), decreasing = TRUE))) +
  theme_bigstatsr(0.9) +
  coord_equal() +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set", color = "True per-variant\nsample size",
       y = "Standard deviations derived from the summary statistics (using maxN)") +
  theme(legend.position = c(0.25, 0.8)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
# ggsave("figures/simu-qc-plot.png", width = 7, height = 7)

# Imputation of N
df_beta3 <- df_beta2
N_est <- (sd_y^2 / sd_val^2 - df_beta2$beta^2) / df_beta2$beta_se^2
N <- df_beta2$n_eff[1]
df_beta3$n_eff <- pmin(pmax(0.5 * N, N_est), 1.1 * N)
plot(df_beta3$n_eff, df_beta$n_eff, pch = 20, col = scales::alpha("black", 0.05))
abline(0, 1, col = "red", lwd = 3)


#### Run lassosum2 and LDpred2-grid with different parameters ####

bigassertr::assert_dir("investigate")

grid <- tidyr::expand_grid(
  method = c("lassosum2", "ldpred2"),
  whichN = c("maxN", "imputeN"),
  corr_num = 1:2
)

grid$res <- purrr::pmap(grid, function(method, whichN, corr_num) {

  cat(".")

  df <- list(maxN = df_beta2, imputeN = df_beta3)[[whichN]]
  corr <- list(corr1, corr2)[[corr_num]]

  res_file <- paste0("investigate/", method, "-mis-",
                     whichN, "-corr", corr_num, ".rds")

  if (method == "lassosum2") {

    runonce::save_run({
      beta_lassosum <- snp_lassosum2(corr, df, ncores = NCORES)
      params2 <- attr(beta_lassosum, "grid_param")
      pred_grid <- big_prodMat(G, beta_lassosum, ind.row = ind.test)
      params2$r2 <- apply(pred_grid, 2, cor, y = y[ind.test])^2
      print(dplyr::arrange(params2, desc(r2)))
      params2
    }, file = res_file)

  } else if (method == "ldpred2") {

    runonce::save_run({
      params <- expand.grid(h2 = signif(0.2 * c(0.01, 0.1, 0.3, 0.7, 1, 1.4), 2),
                            sparse = TRUE, p = signif(seq_log(1e-4, 1, 10), 2))

      beta_grid <- snp_ldpred2_grid(corr, df, params, ncores = NCORES)

      pred_grid <- big_prodMat(G, beta_grid, ind.row = ind.test, ncores = NCORES)
      params$r2 <- apply(pred_grid, 2, cor, y = y[ind.test])^2
      params$sparsity <- colMeans(beta_grid == 0)
      print(dplyr::arrange(params, desc(r2)))
      params
    }, file = res_file)

  } else {
    stop("Method not implemented!")
  }
})


filter(grid, method == "lassosum2") %>%
  mutate(corr = c("normal", "with LD blocks")[corr_num]) %>%
  tidyr::unnest(res) %>%
  ggplot(aes(x = lambda, y = r2, color = paste(corr, whichN, sep = " - "))) +
  facet_wrap(~ delta, labeller = label_bquote(cols = delta: .(delta))) +
  theme_bigstatsr(0.7) +
  geom_point(size = 2) +
  geom_line(size = 1, alpha = 0.5) +
  scale_x_log10(breaks = 10^(-5:0)) +
  scale_color_manual(values = c("#56B4E9", "#999999", "#009E73", "#E69F00")) +
  labs(x = expression(lambda), y = expression(r^2), color = "LD - N") +
  theme(legend.position = "top", legend.key.width = unit(2, "line"))
# ggsave("figures/lassosum2-misN.pdf", width = 9, height = 6)


filter(grid, method == "ldpred2") %>%
  mutate(corr = c("normal", "with LD blocks")[corr_num]) %>%
  tidyr::unnest(res) %>%
  ggplot(aes(x = p, y = r2, color = paste(corr, whichN, sep = " - "))) +
  facet_wrap(~ h2, labeller = label_bquote(cols = h^2: .(h2))) +
  theme_bigstatsr(0.75) +
  geom_point(size = 2) +
  geom_line(size = 1, alpha = 0.5) +
  scale_x_log10(breaks = 10^(-5:0)) +
  scale_color_manual(values = c("#56B4E9", "#999999", "#009E73", "#E69F00")) +
  labs(y = expression(r^2), color = "LD - N") +
  theme(legend.position = "top", legend.key.width = unit(2, "line"),
        panel.spacing = unit(1, "lines"))
# ggsave("figures/ldpred2-misN.pdf", width = 10, height = 6)


#### LDpred2-auto ####

# LDSc reg
(h2_est <- snp_ldsc2(corr0, df_beta)[["h2"]])

corr <- corr2

ldpred2_auto <- snp_ldpred2_auto(corr, df_beta2, h2_init = h2_est,
                                 burn_in = 200, num_iter = 100,
                                 allow_jump_sign = TRUE,
                                 shrink_corr = 1,
                                 verbose = TRUE)

cor(big_prodVec(G, ldpred2_auto[[1]]$beta_est, ind.row = ind.test),
    y[ind.test])^2

plot(ldpred2_auto[[1]]$path_p_est)
plot(ldpred2_auto[[1]]$path_h2_est)

beta_hat <- with(df_beta, beta / sqrt(n_eff * beta_se^2 + beta^2))
plot(ldpred2_auto[[1]]$corr_est, beta_hat); abline(0, 1, col = "red", lwd = 3)
