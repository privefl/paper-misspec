library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- transmute(snp_attach("data/UKBB_HM3_val.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, af_val = freq,
                 a0 = allele1, a1 = allele2)  # reversed somehow..


#### Prostate cancer (PRCA) ####

# download.file("http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip",
#               destfile = "tmp-data/sumstats_PRCA.zip")
# unzip("tmp-data/sumstats_PRCA.zip", exdir = "tmp-data")
sumstats <- fread2(
  "tmp-data/meta_v3_onco_euro_overall_ChrAll_1_release.txt",
  select = c("Chr", "position", "Allele1", "Allele2", "Effect", "StdErr",
             "Pvalue", "Freq1", "OncoArray_imputation_r2"),
  col.names = c("chr", "pos", "a1", "a0", "beta", "beta_se", "p", "freq", "info")
) %>%
  mutate(a0 = toupper(a0), a1 = toupper(a1)) %>%
  filter(info > 0.4, beta_se > 0)

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 13,421,441 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 818,400 variants have been matched; 0 were flipped and 398,171 were reversed.

sizes <- bigreadr::fread2(text = c("44825\t27904","20219\t20440","1854\t1894","3650\t3940",
                                   "474\t482","1458\t512","2068\t2993", "4600\t2941"))
info_snp$n_eff <- print(sum(4 / rowSums(1 / sizes)))  # 135316.1
4 / (1 / 79148 + 1 / 61106)  # 137933.1

mean(info_snp$info)
hist(info_snp$info, "FD")

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
# ggsave("figures/prca_qc.png", width = 11, height = 6)

info_snp2$is_bad2 <- with(info_snp2,
                          sd_ss2 < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
                            sd_ss2 < 0.1 | sd_af < 0.05)

qplot(sd_af, sd_ss2, color = is_bad2, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr(0.8) +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from the allele frequencies",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

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
# Some are reversed?

with(info_snp2, hist(abs(af_val - freq2), "FD", xlim = c(0, 0.2)))

info_snp2$is_bad3 <- with(info_snp2, abs(af_val - freq2) > 0.1)

with(info_snp2, table(is_bad, is_bad2 | is_bad3))
with(info_snp2, table(is_bad2, is_bad3))

# saveRDS(filter(info_snp2, !is_bad),            "data/sumstats/prca_qc1.rds")
# saveRDS(filter(info_snp2, !is_bad2, !is_bad3), "data/sumstats/prca_qc2.rds")
