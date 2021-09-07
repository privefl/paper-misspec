library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- transmute(snp_attach("data/UKBB_HM3_val.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, af_val = freq,
                 a0 = allele1, a1 = allele2)  # reversed somehow..


#### MDD ####

gz <- runonce::download_file(
  "https://pgcdata.med.unc.edu/major_depressive_disorders/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz",
  dir = "tmp-data")
sumstats <- fread2(gz,
                   select = c("CHR", "BP", "A1", "A2", "OR", "SE", "P",
                              "Nca", "Nco", "INFO", "FRQ_A_45396", "FRQ_U_97250", "Neff_half"),
                   col.names = c("chr", "pos", "a1", "a0", "or", "beta_se", "p",
                                 "Nca", "Nco", "info", "freq1", "freq2", "Neff_half")) %>%
  as_tibble() %>%
  mutate(beta = log(or), or = NULL, chr = as.integer(chr),
         n_eff = 2 * Neff_half, Neff_half = NULL,
         n_eff2 = 4 / (1 / Nca + 1 / Nco),
         Nca = NULL, Nco = NULL,
         info = pmin(info, 1)) %>%
  filter(n_eff > (0.5 * max(n_eff)),
         info > 0.4) %>%
  print()

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map))
# 9,343,125 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,049,455 variants have been matched; 0 were flipped and 509,899 were reversed.

hist(info_snp$n_eff, "FD"); abline(v = 4 / (1 / 45396 + 1 / 97250), col = "red")
hist(info_snp$n_eff2, "FD"); abline(v = 4 / (1 / 45396 + 1 / 97250), col = "red")
max(info_snp$n_eff)          # 110463.7
max(info_snp$n_eff2)         # 127996.2
4 / (1 / 45396 + 1 / 97250)  # 123796.3

mean(info_snp$info)
hist(info_snp$info, "FD")

info_snp2 <- info_snp %>%
  mutate(freq = (freq1 * 45396 + freq2 * 97250) / (45396 + 97250),
         sd_af = sqrt(2 * freq * (1 - freq)),
         sd_ss = 2 / sqrt(n_eff * beta_se^2 + beta^2),
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

  qplot(sd_af, sd_ss2, alpha = I(0.5), data = df_sub) +
    theme_bigstatsr(0.7) +
    coord_equal() +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations derived from the allele frequencies",
         y = "Standard deviations derived from the summary statistics\n(divided by sqrt(INFO))"),

  scale = 0.95, nrow = 1, labels = LETTERS[1:2], label_size = 16,
  align = "h", rel_widths = c(1, 1.04)
)
# ggsave("figures/mdd_qc.png", width = 11, height = 6)


diff <- with(info_snp2, sd_af - sd_ss2)
hist(diff, "FD")
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

qplot(af_val, freq2, alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr() +
  coord_equal() +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Allele frequencies from the validation set",
       y = "Allele frequencies from the summary statistics")
# Some are reversed?

with(info_snp2, hist(abs(af_val - freq2), "FD", xlim = c(0, 0.2)))

info_snp2$is_bad3 <- with(info_snp2, abs(af_val - freq2) > 0.1)

with(info_snp2, table(is_bad, is_bad2 | is_bad3))
with(info_snp2, table(is_bad2, is_bad3))

# saveRDS(filter(info_snp2, !is_bad),            "data/sumstats/mdd_qc1.rds")
# saveRDS(filter(info_snp2, !is_bad2, !is_bad3), "data/sumstats/mdd_qc2.rds")
