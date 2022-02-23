library(bigsnpr)
library(ggplot2)
theme_set(theme_bigstatsr())

simu_true <- snp_attach("data/ukbb4simu.rds")
G_true <- simu_true$genotypes
INFO0 <- simu_true$map$info   # the ones reported in UKBB
NCORES <- nb_cores()

pheno <- snp_simuPheno(G_true, h2 = 0.2, M = 2000, ncores = NCORES)
y <- pheno$pheno

load("data/ukbb4simu_ind.RData")
gwas_true <- big_univLinReg(G_true, y[ind.gwas], ind.train = ind.gwas, ncores = NCORES)

simu_imp <- snp_attach("data/ukbb4simu_imp.rds")
INFO <- simu_imp$map$info  # the ones recomputed for the subset
qplot(INFO, INFO0) +
  geom_abline(col = "red", lwd = 2) +
  coord_equal() +
  labs(x = "INFO scores recomputed from the European subset",
       y = "INFO scores reported in the UK Biobank")
# ggsave("figures/compare-info.png", width = 7, height = 7)

G_imp <- simu_imp$genotypes
gwas_imp <- big_univLinReg(G_imp, y[ind.gwas], ind.train = ind.gwas, ncores = NCORES)

scale_true <- big_scale()(G_true, ind.row = ind.gwas)
scale_imp <- big_scale()(G_imp, ind.row = ind.gwas)
# scale_freq <- sqrt(scale_imp$center * (1 - scale_imp$center / 2))

source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
plot_grid2(list(
  qplot(scale_true$scale, scale_imp$scale, color = INFO) +
    geom_abline(color = "red") +
    coord_equal() +
    scale_color_viridis_c() +
    guides(color = guide_colorbar(barheight = 15, ticks.linewidth = 2)) +
    labs(x = "Standard deviations of \"true\" genotypes",
         y = "Standard deviations of imputed dosages"),
  qplot(scale_true$scale, scale_imp$scale / sqrt(INFO), color = INFO) +
    geom_abline(color = "red") +
    coord_equal() +
    scale_color_viridis_c() +
    labs(x = "Standard deviations of \"true\" genotypes",
         y = "Corrected standard deviations of imputed dosages")
), scale = 0.95, nrow = 1, labels = c("A", "B"), label_size = 18,
title_ratio = 0, legend_ratio = 0.1)
# ggsave("figures/compare-sd.png", width = 12.5, height = 5.5, scale = 1.3)

plot_grid2(list(
  qplot(gwas_imp$estim, gwas_true$estim, color = INFO) +
    geom_abline(color = "red") +
    scale_color_viridis_c() +
    guides(color = guide_colorbar(barheight = 15, ticks.linewidth = 2)) +
    labs(x = "GWAS effect sizes from imputed dosages",
         y = "GWAS effect sizes from \"true\" genotypes"),
  qplot(gwas_imp$estim * sqrt(INFO), gwas_true$estim, color = INFO) +
    geom_abline(color = "red") +
    scale_color_viridis_c() +
    labs(x = "Corrected GWAS effect sizes from imputed dosages",
         y = "GWAS effect sizes from \"true\" genotypes")
), scale = 0.95, nrow = 1, labels = c("A", "B"), label_size = 18,
title_ratio = 0, legend_ratio = 0.1)
# ggsave("figures/compare-beta.png", width = 12.5, height = 5.5, scale = 1.3)

plot_grid2(list(
  qplot(gwas_imp$std.err, gwas_true$std.err, color = INFO) +
    geom_abline(color = "red") +
    scale_color_viridis_c() +
    guides(color = guide_colorbar(barheight = 15, ticks.linewidth = 2)) +
    labs(x = "SEs of GWAS effect sizes from imputed dosages",
         y = "SEs of GWAS effect sizes from \"true\" genotypes"),
  qplot(gwas_imp$std.err * sqrt(INFO), gwas_true$std.err, color = INFO) +
    geom_abline(color = "red") +
    scale_color_viridis_c() +
    labs(x = "Corrected SEs of GWAS effect sizes from imputed dosages",
         y = "SEs of GWAS effect sizes from \"true\" genotypes")
), scale = 0.95, nrow = 1, labels = c("A", "B"), label_size = 18,
title_ratio = 0, legend_ratio = 0.1)
# ggsave("figures/compare-beta-se.png", width = 12.5, height = 5.5, scale = 1.3)


#### Multiple imputation ####

ind <- seq(1, ncol(G_true), by = 10)
snp_id <- with(simu_true$map[ind, ],
               paste(chromosome, physical.pos, allele1, allele2, sep = "_"))
d <- 20

system.time({
  tmp <- tempfile(tmpdir = "tmp-data")
  rds <- snp_readBGEN(
    bgenfiles = "UKBB/bgen/ukb_imp_chr22_v3.bgen",
    list_snp_id = list(rep(snp_id, each = d)),
    backingfile = tmp,
    ind_row = simu_true$fam$ind_bgen[ind.gwas],
    ncores = NCORES,
    read_as = "random"
  )
  G <- snp_attach(rds)$genotypes
  gwas_samp <- big_univLinReg(G, y[ind.gwas], ncores = NCORES)
  file.remove(paste0(tmp, c(".bk", ".rds")))
})
# by 40 -> 5 min
# by 10 -> 20 min

#  https://doi.org/10.1371/journal.pgen.1006091
all_betas <- matrix(gwas_samp$estim, ncol = d, byrow = TRUE)
beta_MI <- rowMeans(all_betas)
s2_W <- rowMeans(matrix(gwas_samp$std.err^2, ncol = d, byrow = TRUE))
s2_B <- rowSums(sweep(all_betas, 1, beta_MI, '-')^2) / (d - 1)
s2_MI <- s2_W + (1 + 1 / d) * s2_B
beta_se_MI <- sqrt(s2_MI)
df_MI <- (d - 1) * (1 + d * s2_W / ((d + 1) * s2_B))^2
lpval_MI <- -(log(2) + pt(abs(beta_MI / beta_se_MI), df = df_MI,
                              lower.tail = FALSE, log.p = TRUE)) / log(10)
hist(10^-lpval_MI)
lpval_true <- -predict(gwas_true)
qplot(lpval_MI, lpval_true[ind], color = INFO[ind], alpha = I(0.8)) +
  geom_abline(color = "red") +
  scale_color_viridis_c()

SEQ <- seq(0.1, 0.5, by = 0.01)
diff <- sapply(SEQ, function(power) crossprod(beta_se_MI - gwas_imp$std.err[ind] * INFO[ind]^power))
plot(SEQ, diff)
qplot(beta_se_MI, gwas_imp$std.err[ind] * INFO[ind]^0.17, color = INFO[ind]) +
  geom_abline(col = "red") +
  scale_color_viridis_c() +
  labs(x = "beta_se_MI", y = "beta_se_imp", color = "INFO")
qplot(beta_MI / beta_se_MI, (gwas_true$estim / gwas_true$std.err)[ind], color = INFO[ind]) +
  geom_abline(col = "red") +
  scale_color_viridis_c() +
  labs(x = "Z_MI", y = "Z_true", color = "INFO")
qplot(beta_MI / beta_se_MI, (gwas_true$estim / gwas_true$std.err * sqrt(INFO))[ind], color = INFO[ind]) +
  geom_abline(col = "red") +
  scale_color_viridis_c() +
  labs(x = "Z_MI", y = "Z_true * sqrt(INFO)", color = "INFO")

plot_grid2(list(
  qplot(beta_MI / beta_se_MI, (gwas_imp$estim / gwas_imp$std.err * (INFO))[ind], color = INFO[ind]) +
    geom_abline(col = "red") +
    scale_color_viridis_c() +
    guides(color = guide_colorbar(barheight = 15, ticks.linewidth = 2)) +
    labs(x = "Z-scores from multiple imputation",
         y = "Overcorrected Z-scores from imputed dosages",
         color = "INFO"),
  qplot(beta_MI, (gwas_imp$estim * INFO)[ind], color = INFO[ind]) +
    geom_abline(color = "red") +
    scale_color_viridis_c() +
    labs(x = "GWAS effect sizes from multiple imputation",
         y = "Overcorrected GWAS effect sizes from imputed dosages")
), scale = 0.95, nrow = 1, labels = c("A", "B"), label_size = 18,
title_ratio = 0, legend_ratio = 0.1)
# ggsave("figures/compare-mult-imp.png", width = 12.5, height = 5.5, scale = 1.3)
