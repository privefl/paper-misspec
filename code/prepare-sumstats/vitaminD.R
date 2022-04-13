library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)


map <- transmute(snp_attach("data/UKBB_HM3_val.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, rsid, af_val = freq,
                 a0 = allele1, a1 = allele2)  # reversed somehow..


# See google drive link in https://doi.org/10.1038/s41467-017-02662-2
# Should have used Allele1 -> a1 and Allele2 -> a0
sumstats <- fread2("tmp-data/25HydroxyVitaminD_QC.METAANALYSIS1.txt",
                   select = c("Chr", "Pos", "MarkerName", "Allele1", "Allele2",
                              "Effect", "StdErr", "SAMPLESIZE", "HetPVal"),
                   col.names = c("chr", "pos", "rsid", "a0", "a1", "beta", "beta_se", "N", "het_pval")) %>%
  mutate_at(c("a0", "a1"), toupper)
head(sumstats$pos)
# 127246000  33550000  24178800  96896300 160770000  59208100  -> useless for matching..


info_snp <- as_tibble(snp_match(sumstats, map, join_by_pos = FALSE))
# 2,579,296 variants to be matched.
# 2 ambiguous SNPs have been removed.
# 1,016,935 variants have been matched; 6,394 were flipped and 522,741 were reversed.

hist(-log10(info_snp$het_pval), "FD")

info_snp2 <- info_snp %>%
  filter(het_pval > 1e-5) %>%
  mutate(n_eff = N,
         sd_af = sqrt(2 * af_val * (1 - af_val)),
         sd_ss = 1 / sqrt(n_eff * beta_se^2 + beta^2),
         sd_ss = sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5))


info_snp2$is_bad <- with(info_snp2,
                         sd_ss < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
                           sd_ss < 0.1 | sd_af < 0.05)

qplot(sd_af, sd_ss, color = ifelse(is_bad, "Yes", "No"), alpha = I(0.5),
      data = slice_sample(info_snp2, n = 100e3)) +
  theme_bigstatsr(0.8) +
  coord_equal() +
  theme(legend.position = c(0.2, 0.8)) +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")
# ggsave("figures/vitaminD_qc.png", width = 7.5, height = 8.5)


# pdf("figures/hist_N_vitaminD.pdf", width = 7, height = 4.5)
hist(info_snp$N, breaks = 50, main = NULL, xlab = "Sample size")
abline(v = 0.7 * max(info_snp$N), col = "red")
# dev.off()

# saveRDS(info_snp2, "data/sumstats/vitaminD_noqc.rds")
# saveRDS(filter(info_snp2, !is_bad), "data/sumstats/vitaminD_qc1.rds")
# saveRDS(filter(info_snp2, !is_bad, N > (0.7 * max(N))), "data/sumstats/vitaminD_qc2.rds")
