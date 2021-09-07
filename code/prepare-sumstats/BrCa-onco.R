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

# OncoArray
sumstats <- fread2("tmp-data/sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_onco2_eaf_controls",
                              "bcac_onco2_beta",
                              "bcac_onco2_se",
                              "bcac_onco2_r2"),
                   col.names = c("chr", "pos", "a0", "a1", "freq",
                                 "beta", "beta_se", "info")) %>%
  filter(info > 0.4, beta_se > 0)

Neff1 <- 4 / (1 / 61282 + 1 / 45494)  # Main text

# Supp Table 1
Ncases <- c(1068, 589, 143, 848, 0, 676, 100, 782, 454, 137,
            102, 1437, 1378, 265, 678, 681, 11, 1406, 3029, 1100, 646, 3535,
            3, 3365, 459, 385, 357, 916, 426, 53, 0,
            0, 70, 221, 77, 2205, 136, 0, 805, 90, 57, 788, 651, 615, 70,
            692, 375, 44, 0, 1283, 641, 0, 776, 1588, 1606, 993, 1015, 1439, 865,
            1019, 968, 451, 0, 3999, 0, 0, 2014, 1017, 1478,
            439, 1335, 387, 759, 0, 501, 1159, 0, 4908)
sum(Ncases) # 62565 -> 1283 more than in the text (because NBCS study was excluded)
Ncontrols <- c(0, 188, 4, 374, 0, 248, 442, 834, 27, 0,
               0, 723, 725, 167, 817, 332, 3, 712, 3025, 577, 0, 3522,
               3, 1593, 284, 0, 180, 865, 0, 2, 0,
               0, 214, 131, 0, 6042, 182, 0, 435, 93, 288, 366, 179, 712, 127,
               1523, 1605, 29, 0, 0, 613, 0, 149, 1804, 1905, 217, 660, 1658, 858,
               0, 0, 231, 0, 989, 0, 0, 1560, 0, 708, 0, 0, 157, 0, 0, 258, 567, 974, 4613)
sum(Ncontrols) # 45494

sumstats$n_eff <- Neff1

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 11,726,603 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,054,233 variants have been matched; 0 were flipped and 0 were reversed.

# pdf("figures/brca_onco_hist_info.pdf", width = 7, height = 6)
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
    theme(legend.position = c(0.18, 0.68)),

  qplot(sd_af, sd_ss2, alpha = I(0.5), data = df_sub,
        color = ifelse(chr %in% c(6, 8), "Chromosomes 6 and 8", "Other chromosomes")) +
    theme_bigstatsr(0.7) +
    coord_equal() +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations derived from the allele frequencies",
         y = "Standard deviations derived from the summary statistics\n(divided by sqrt(INFO))",
         color = NULL) +
    scale_color_manual(values = c("#000000", "#E69F00")) +
    theme(legend.position = c(0.28, 0.85)),

  scale = 0.95, nrow = 1, labels = LETTERS[1:2], label_size = 16,
  align = "h", rel_widths = c(1, 1.04)
)
# ggsave("figures/brca_onco_qc.png", width = 11, height = 6)

diff <- with(info_snp2, sd_af - sd_ss2)
hist(abs(diff), "FD", xlim = c(0, 0.1))
diff2 <- bigutilsr::rollmean(diff, 10)
hist(diff2, "FD", xlim = c(0, 0.1))

qplot(pos / 1e6, data = info_snp2[diff2 > 0.03, ], breaks = 5:35,
      fill = as.factor(chr), color = I("#000000"), alpha = I(0.6)) +
  scale_x_continuous(breaks = seq(0, 1000, by = 5)) +
  theme_bigstatsr(0.75) +
  labs(fill = "Chr", x = "Physical position (in Mbp)", y = "Frequency") +
  theme(legend.position = c(0.5, 0.75)) +
  scale_fill_manual(values = c("#009E73", "#CC79A7"))
# ggsave("figures/hist_bad_pos.pdf", width = 7, height = 5)

info_snp2$is_bad2 <- (diff2 > 0.03 | info_snp2$sd_ss2 < 0.1 | info_snp2$sd_af < 0.05)


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

with(info_snp2, hist(abs(af_val - freq2), "FD", xlim = c(0, 0.2)))

sum(info_snp2$is_bad3 <- with(info_snp2, abs(af_val - freq2) > 0.1))  # 21

with(info_snp2, table(is_bad2, is_bad3))  # none in common
with(info_snp2, table(is_bad, is_bad2 | is_bad3))

# saveRDS(filter(info_snp2, !is_bad),            "data/sumstats/brca_onco_qc1.rds")
# saveRDS(filter(info_snp2, !is_bad2, !is_bad3), "data/sumstats/brca_onco_qc2.rds")
