library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- transmute(snp_attach("data/UKBB_HM3_val.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, af_val = freq,
                 a0 = allele1, a1 = allele2)  # reversed somehow..


#### Breast cancer (BRCA) ####

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "tmp-data/sumstats_BRCA.txt.gz")
# R.utils::gunzip("tmp-data/sumstats_BRCA.txt.gz")

# iCOGS
sumstats <- fread2("tmp-data/sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_icogs2_eaf_controls",
                              "bcac_icogs2_beta",
                              "bcac_icogs2_se",
                              "bcac_icogs2_r2"),
                   col.names = c("chr", "pos", "a0", "a1", "freq",
                                 "beta", "beta_se", "info")) %>%
  filter(info > 0.4, beta_se > 0)

Neff1 <- 4 / (1 / 46785  + 1 / 42892)
Neff2 <- 4 / (1 / 45290  + 1 / 41880)
Neff1 / Neff2  # 1.03
sumstats$n_eff <- Neff2

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 11,417,746 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,051,242 variants have been matched; 0 were flipped and 0 were reversed.

# pdf("figures/brca_icogs_hist_info.pdf", width = 7, height = 6)
# par(mar = c(4.1, 4.1, 0.1, 0.1))
# hist(info_snp$info, main = NULL, xlab = "INFO")
# dev.off()
mean(info_snp$info)

info_snp2 <- mutate(info_snp,
                    sd_af = sqrt(2 * freq * (1 - freq)),
                    sd_ss = 2 / sqrt(n_eff * beta_se^2 + beta^2),
                    sd_ss2 = sd_ss / sqrt(info))

# original filtering applied in LDpred2 paper (but using af from sumstats)
info_snp2$is_bad <- with(info_snp2,
                         sd_ss < (0.5 * sd_af) | sd_ss > (sd_af + 0.1) |
                           sd_ss < 0.1 | sd_af < 0.05)
mean(.Last.value)

qplot(sd_af, sd_ss, color = is_bad, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr(0.8) +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from the allele frequencies",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

# impact of INFO score
df_sub <- slice_sample(info_snp2, n = 100e3)
plot_grid(
  qplot(sd_af, sd_ss, color = info, alpha = I(0.5), data = df_sub) +
    theme_bigstatsr(0.7) +
    coord_equal() +
    scale_color_viridis_c(direction = -1) +
    guides(color = guide_colorbar(barheight = 10, ticks.linewidth = 1.5)) +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations derived from the allele frequencies",
         y = "Standard deviations derived from the summary statistics",
         color = "INFO") +
    theme(legend.position = c(0.2, 0.7)),

  qplot(sd_af, sd_ss2, alpha = I(0.5), data = df_sub) +
    theme_bigstatsr(0.7) +
    coord_equal() +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations derived from the allele frequencies",
         y = "Standard deviations derived from the summary statistics\n(divided by sqrt(INFO))"),

  scale = 0.95, nrow = 1, labels = LETTERS[1:2], label_size = 16,
  align = "h", rel_widths = c(1, 1.04)
)
# ggsave("figures/brca_icogs_qc.png", width = 11, height = 6)

diff <- with(info_snp2, sd_af - sd_ss2)
hist(abs(diff), "FD")

info_snp2$is_bad2 <- (info_snp2$sd_ss2 < 0.1 | info_snp2$sd_af < 0.05)


# comparing allele frequencies
info_snp2$freq2 <- ifelse(info_snp2$beta * sumstats$beta[info_snp2$`_NUM_ID_.ss`] < 0,
                          1 - info_snp2$freq, info_snp2$freq)

qplot(af_val, freq2, color = is_bad2, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

with(info_snp2, hist(af_val - freq2, "FD", xlim = c(-0.2, 0.2)))

sum(info_snp2$is_bad3 <- with(info_snp2, abs(af_val - freq2) > 0.1))  # 459

with(info_snp2, table(is_bad2, is_bad3))
with(info_snp2, table(is_bad, is_bad2 | is_bad3))

# saveRDS(filter(info_snp2, !is_bad),            "data/sumstats/brca_icogs_qc1.rds")
# saveRDS(filter(info_snp2, !is_bad2, !is_bad3), "data/sumstats/brca_icogs_qc2.rds")
