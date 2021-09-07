library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- transmute(snp_attach("data/UKBB_HM3_val.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, af_val = freq,
                 a0 = allele1, a1 = allele2)  # reversed somehow..


#### Type 1 diabetes (T1D) ####

# To request a download link, go to
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.ns8q3
# download.file("http://merritt.cdlib.org/cloudcontainer/mrtstore2/35227889.tar.gz",
#               destfile = "tmp-data/sumstats_T1D.tar.gz")
# untar("tmp-data/sumstats_T1D.tar.gz", exdir = "T1D")

# Affymetrix
sumstats <- fread2(
  paste0("../paper-ldpred2/tmp-data/meta_chr_", 1:22),
  select = c("chromosome", "position", "a0", "a1", "info_score.A",
             "EUR_MAF_1kG", "controls_maf.A", "aff.OR", "aff.se", "qc.check"),
  col.names = c("chr", "pos", "a0", "a1", "info",
                "freq_1000G", "freq", "OR", "beta_se", "qc")) %>% # not freq, but MAF
  mutate(beta = log(OR), OR = NULL) %>%
  filter(info > 0.4, beta_se > 0)

# https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002362#pmed-1002362-t001
sumstats$n_eff <- 4 / (1 / 1930 + 1 / (1455 + 1490 + 1884))

(info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map)))
# 8,552,620 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 934,712 variants have been matched; 520 were flipped and 1,763 were reversed.

table(info_snp$qc)
mean(info_snp$info)
hist(info_snp$info, "FD")
summary(info_snp$info)

info_snp2 <- mutate(info_snp,
                    sd_val = sqrt(2 * af_val * (1 - af_val)),
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

  qplot(sd_af, sd_ss2, color = as.factor(qc), alpha = I(0.5), data = df_sub) +
    theme_bigstatsr(0.7) +
    coord_equal() +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations derived from the allele frequencies",
         y = "Standard deviations derived from the summary statistics\n(divided by sqrt(INFO))",
         color = "QC") +
    theme(legend.position = c(0.25, 0.75)),

  scale = 0.95, nrow = 1, labels = LETTERS[1:2], label_size = 16,
  align = "h", rel_widths = c(1, 1.04)
)
# ggsave("figures/t1d_affy_qc.png", width = 11, height = 6)

diff <- with(info_snp2, sd_af - sd_ss2)
hist(abs(diff), "FD", xlim = c(0, 0.1))
info_snp2$is_bad2 <- (abs(diff) > 0.1 | info_snp2$sd_ss2 < 0.1 | info_snp2$sd_af < 0.05)

qplot(sd_af, sd_ss2, color = is_bad2, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr() +
  coord_equal() +
  geom_abline(linetype = 2, color = "red") +
  scale_color_viridis_d(direction = -1) +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

# original filtering applied in LDpred2 paper (using af from validation)
diff2 <- with(info_snp2, sd_val - sd_ss2)
hist(abs(diff2), "FD", xlim = c(0, 0.1))
info_snp2$is_bad3 <- (abs(diff2) > 0.1 | info_snp2$sd_val < 0.05)

qplot(sd_val, sd_ss2, color = is_bad3, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr(0.8) +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from the allele frequencies",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

with(info_snp2, table(is_bad2, is_bad3))
with(info_snp2, table(is_bad, is_bad2 | is_bad3))

# saveRDS(filter(info_snp2, !is_bad),            "data/sumstats/t1d_affy_qc1.rds")
# saveRDS(filter(info_snp2, !is_bad2, !is_bad3), "data/sumstats/t1d_affy_qc2.rds")
