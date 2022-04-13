library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- transmute(snp_attach("data/UKBB_HM3_val.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, rsid,
                 a0 = allele1, a1 = allele2)  # reversed somehow..

zip <- runonce::download_file(
  "https://pheweb.jp/download/Height",
  dir = "tmp-data", fname = "sumstats_height_BBJ.zip")
unzip(zip, exdir = "tmp-data", overwrite = FALSE)

sumstats <- bigreadr::fread2(
  "tmp-data/hum0197.v3.BBJ.Hei.v1/GWASsummary_Height_Japanese_SakaueKanai2020.auto.txt.gz",
  select = c("CHR", "BP", "SNP", "ALLELE0", "ALLELE1", "A1FREQ", "INFO",
             "BETA", "SE", "CHISQ_LINREG", "CHISQ_BOLT_LMM_INF"),
  col.names = c("chr", "pos", "rsid", "a0", "a1", "freq", "info",
                "beta", "beta_se", "CHISQ_LINREG", "CHISQ_BOLT_LMM_INF")
) %>%
  filter(chr %in% 1:22, pmin(freq, 1 - freq) > 0.01) %>%
  mutate(chr = as.integer(chr)) %>%
  print()

# Estimate power gain using the random effects
qplot(CHISQ_BOLT_LMM_INF / CHISQ_LINREG, color = as.factor(chr), geom = "density",
      data = filter(sumstats, CHISQ_LINREG > 30))

power <-
  with(filter(sumstats, CHISQ_LINREG > 30),
       matrixStats::weightedMedian(CHISQ_BOLT_LMM_INF / CHISQ_LINREG,
                                   w = CHISQ_BOLT_LMM_INF))
# 1.143

info_snp <- as_tibble(bigsnpr::snp_match(sumstats, map, join_by_pos = FALSE))
# 7,444,297 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 870,000 variants have been matched; 0 were flipped and 870,000 were reversed.

hist(info_snp$info)
mean(info_snp$info)  # 98.6%

info_snp2 <- info_snp %>%
  mutate(n_eff = 165056 * power,
         sd_af = sqrt(2 * freq * (1 - freq)),
         sd_ss = 1 / sqrt(n_eff * beta_se^2 + beta^2),
         sd_ss2 = sd_ss / sqrt(info),
         sd_ss2 = sd_ss2 / quantile(sd_ss2, 0.99) * sqrt(0.5))

qplot(sd_af, sd_ss, color = info, alpha = I(0.2),
      data = slice_sample(info_snp2, n = 50e3)) +
  theme_bw(15) +
  coord_equal() +
  scale_color_viridis_c(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics")

qplot(sd_af, sd_ss2, color = info, alpha = I(0.2),
      data = slice_sample(info_snp2, n = 50e3)) +
  theme_bw(15) +
  coord_equal() +
  scale_color_viridis_c(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics")


with(info_snp2, hist(sd_af - sd_ss2))

saveRDS(filter(info_snp2, abs(sd_af - sd_ss2) < 0.03),
        "data/sumstats_BBJ/height.rds")
