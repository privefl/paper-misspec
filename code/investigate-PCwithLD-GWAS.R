library(bigsnpr)
ukb <- snp_attach("data/UKBB_HM3_val.rds")
G <- ukb$genotypes
PC19 <- as.matrix(bigreadr::fread2("UKBB/ukb41181.csv", select = "22009-0.19")[ukb$fam$id_csv, 1])

ind_chr6 <- which(ukb$map$chromosome == 6)
y <- snp_simuPheno(G, h2 = 0.2, M = 2000, ind.possible = ind_chr6)$pheno

POS <- ukb$map$physical.pos[ind_chr6]
is_in_bad_region <- (POS > 70e6) & (POS < 91e6)

gwas  <- big_univLinReg(G, y, ind.col = ind_chr6, ncores = nb_cores())
gwas2 <- big_univLinReg(G, y, ind.col = ind_chr6, ncores = nb_cores(), covar.train = PC19)

sd_ss  <- with(gwas,  1 / sqrt(length(y) * std.err^2 + estim^2))
sd_ss2 <- with(gwas2, 1 / sqrt(length(y) * std.err^2 + estim^2))
sd_af <- with(ukb$map[ind_chr6, ], sqrt(2 * freq * (1 - freq)))
INFO <- ukb$map$info[ind_chr6]

library(ggplot2)
plot_grid(
  qplot(sd_af, sd_ss / sqrt(INFO), color = is_in_bad_region, alpha = I(0.5)) +
    theme_bigstatsr(0.7) +
    coord_equal() +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations derived from the allele frequencies",
         y = "Standard deviations derived from the summary statistics",
         color = "Is in 70-91 Mbp\nof chromosome 6?") +
    scale_color_manual(values = c("#E69F00", "#000000")) +
    theme(legend.position = c(0.28, 0.8)),

  qplot(sd_af, sd_ss2 / sqrt(INFO), color = is_in_bad_region, alpha = I(0.5)) +
    theme_bigstatsr(0.7) +
    coord_equal() +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations derived from the allele frequencies",
         y = "Standard deviations derived from the summary statistics",
         color = "Is in 70-91 Mbp\nof chromosome 6?") +
    scale_color_manual(values = c("#E69F00", "#000000")) +
    theme(legend.position = c(0.28, 0.8)),

  scale = 0.95, nrow = 1, labels = LETTERS[1:2], label_size = 16
)
# ggsave("figures/gwas_bad_pc.png", width = 12, height = 6)
