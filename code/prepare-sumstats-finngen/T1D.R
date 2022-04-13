library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- transmute(snp_attach("data/UKBB_HM3_val.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, rsid,
                 a0 = allele1, a1 = allele2)  # reversed somehow..

# https://www.finngen.fi/en/access_results
file <- R.utils::gunzip("tmp-data/summary_stats_finngen_R6_E4_DM1.gz",
                        skip = TRUE, remove = FALSE)

writeLines(readLines(file, n = 2))

sumstats <- bigreadr::fread2(
  file,
  select = c("#chrom", "pos", "rsids", "ref", "alt", "af_alt", "beta", "sebeta",
             "n_hom_cases", "n_hom_ref_cases", "n_het_cases",
             "n_hom_controls", "n_hom_ref_controls", "n_het_controls"),
  col.names = c("chr", "pos", "rsid", "a0", "a1", "freq", "beta", "beta_se",
                "n_hom_cases", "n_hom_ref_cases", "n_het_cases",
                "n_hom_controls", "n_hom_ref_controls", "n_het_controls"))

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map, join_by_pos = FALSE))
# 16,355,179 variants to be matched.
# 246 ambiguous SNPs have been removed.
# 1,042,134 variants have been matched; 690 were flipped and 2,804 were reversed.

info_snp2 <- info_snp %>%
  filter(chr %in% 1:22, pmin(freq, 1 - freq) > 0.01) %>%
  mutate(n_case = round(n_hom_cases + n_hom_ref_cases + n_het_cases),
         n_control = round(n_hom_controls + n_hom_ref_controls + n_het_controls))

table(info_snp2$n_case)    # all 7609
table(info_snp2$n_control) # all 215,160

info_snp3 <- info_snp2 %>%
  select(!starts_with("n_")) %>%
  mutate(n_eff = 4 / (1 / 7609 + 1 / 215160),
         sd_af = sqrt(2 * freq * (1 - freq)),
         sd_ss = 2 / sqrt(n_eff * beta_se^2 + beta^2))

qplot(sd_af, sd_ss, alpha = I(0.2),
      data = slice_sample(info_snp3, n = 50e3)) +
  theme_bw(15) +
  coord_equal() +
  scale_color_viridis_c(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics")

n_eff_imp <- with(info_snp3, (4 / sd_af^2 - beta^2) / beta_se^2)
Neff <- info_snp3$n_eff[1]; hist(n_eff_imp, "FD", xlim = c(0, Neff)); abline(v = Neff, col = "red")
new_Neff <- quantile(n_eff_imp, 0.8)
new_Neff / Neff # 73.1%

info_snp4 <- info_snp3 %>%
  mutate(n_eff = new_Neff,
         sd_ss = 2 / sqrt(n_eff * beta_se^2 + beta^2))

qplot(sd_af, sd_ss, alpha = I(0.2),
      data = slice_sample(info_snp4, n = 50e3)) +
  theme_bw(15) +
  coord_equal() +
  scale_color_viridis_c(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics")

with(info_snp4, hist(sd_af - sd_ss))

saveRDS(filter(info_snp4, abs(sd_af - sd_ss) < 0.03),
        "data/sumstats_FIN/t1d.rds")
