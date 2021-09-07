library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- transmute(snp_attach("data/UKBB_HM3_val.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, af_val = freq,
                 a0 = allele1, a1 = allele2)  # reversed somehow..


#### Coronary artery disease (CAD) ####

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
#               destfile = "tmp-data/sumstats_CAD.txt")
# bigreadr::nlines("tmp-data/sumstats_CAD.txt")  # 9,455,779
sumstats <- fread2("tmp-data/sumstats_CAD.txt",
                   select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "se_dgc", "p_dgc",
                              "effect_allele_freq", "median_info", "n_studies", "het_pvalue"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se",
                                 "p", "freq", "info", "n_studies", "het_pvalue"))

study_size <- c(
  # from Supp Table 1
  "5719 6545","206 259","278 312","505 1021","392 410","1010 3998",
  "1628 368","2083 2048","1216 653","658 5841","1802 466","2099 2690",
  "634 1608","1207 1288","1061 1467","1089 1147","877 2187",
  "2700 2758","361 2778","487 1381","758 3337","2791 3757","2095 503",
  "933 468","2905 2998","947 1008","1294 1529","843 318","933 468",
  "119 830","631 334","836 761","426 594","814 5999","322 857",
  "1926 2938","4651 4452","4380 3929","1535 772","1007 22286","402 448",
  "745 1389","397 2474","506 5335","259 4202","334 3446","2034 3210","454 8443")
study_size2 <- fread2(text = study_size, col.names = c("Nca", "Nco"))
nrow(study_size2)   # 48
colSums(study_size2)   # Nca: 61289 -- Nco: 126310  (slightly different from 60801 -- 123504)
(Neff <- sum(4 / rowSums(1 / study_size2)))  # 129014.3
4 / (1 / 60801 + 1 / 123504)  # 162972.6


info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 9,455,778 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,052,200 variants have been matched; 0 were flipped and 730,526 were reversed.

mean(info_snp$info)
hist(info_snp$info, "FD")
hist(info_snp$n_studies); abline(v = 48, col = "red")
hist(-log10(info_snp$het_pvalue), "FD")


info_snp2 <- info_snp %>%
  filter(n_studies > 35) %>%
  mutate(sd_af = sqrt(2 * freq * (1 - freq)),
         sd_ss = 2 / sqrt(Neff * beta_se^2 + beta^2),
         sd_ss2 = sd_ss / sqrt(info))


info_snp2$is_bad <- with(info_snp2,
                         sd_ss < (0.5 * sd_af) | sd_ss > (sd_af + 0.1) |
                           sd_ss < 0.1 | sd_af < 0.05)

qplot(sd_af, sd_ss, color = is_bad, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
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

  qplot(sd_af, sd_ss2, color = n_studies, alpha = I(0.5), data = df_sub) +
    theme_bigstatsr(0.7) +
    coord_equal() +
    geom_abline(linetype = 2, color = "red") +
    scale_color_viridis_c(direction = -1) +
    labs(x = "Standard deviations derived from the allele frequencies",
         y = "Standard deviations derived from the summary statistics\n(divided by sqrt(INFO))",
         color = "#studies") +
    theme(legend.position = c(0.25, 0.75)),

  scale = 0.95, nrow = 1, labels = LETTERS[1:2], label_size = 16,
  align = "h", rel_widths = c(1, 1.04)
)
# ggsave("figures/cad_qc.png", width = 11, height = 6)


info_snp2$is_bad2 <- with(info_snp2,
                          sd_ss2 < (0.7 * sd_af) | sd_ss2 > (sd_af + 0.1) |
                            sd_ss2 < 0.1 | sd_af < 0.05)

qplot(sd_af, sd_ss2, color = is_bad2, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

qplot(sd_af, sd_ss2, color = het_pvalue < 1e-7, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from validation set",
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

with(info_snp2, hist(abs(af_val - freq2), "FD", xlim = c(0, 0.2)))

info_snp2$is_bad3 <- with(info_snp2, abs(af_val - freq2) > 0.15)

with(info_snp2, table(is_bad2, het_pvalue < 1e-7))  # all het captured
with(info_snp2, table(is_bad2, is_bad3))            # almost none in common

info_snp2$n_eff <- with(info_snp2, Neff * n_studies / max(n_studies))

# saveRDS(filter(info_snp2, !is_bad),            "data/sumstats/cad_qc1.rds")
# saveRDS(filter(info_snp2, !is_bad2, !is_bad3), "data/sumstats/cad_qc2.rds")
