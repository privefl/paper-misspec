file <- runonce::download_file(
  "ftp://share.sph.umich.edu/UKBB_SAIGE_HRC/PheCode_191.11_SAIGE_MACge20.txt.vcf.gz",
  dir = "tmp-data"
)
R.utils::gunzip(file)
sumstats <- bigreadr::fread2("tmp-data/PheCode_191.11_SAIGE_MACge20.txt.vcf")
sumstats <- dplyr::filter(sumstats, af > 0.005)

hm3_info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))

match_info <- bigsnpr::snp_match(
  sumstats = dplyr::rename(sumstats, chr = "#CHROM", pos = POS, a1 = ALT, a0 = REF),
  info_snp = hm3_info
)
# 11,904,913 variants to be matched.
# 1,790,164 ambiguous SNPs have been removed.
# 1,119,564 variants have been matched; 0 were flipped and 783,941 were reversed.

# https://github.com/privefl/paper-ldpred2/blob/master/paper/paper-ldpred2-supp.pdf
df <- dplyr::mutate(match_info,
                    sd_af = sqrt(2 * af * (1 - af)),
                    n_eff = 4 / (1 / num_cases + 1 / num_controls),
                    sd_ss = 2 / (sebeta * sqrt(n_eff)))

# Reduce points for plotting (but try to keep outliers)
hist(ratio <- with(df, sd_ss / sd_af))
hist(out <- abs(ratio - median(ratio)) / mad(ratio))
ind <- sample(nrow(df), 50e3, prob = out)

library(ggplot2)
theme_set(bigstatsr::theme_bigstatsr())

ggplot(df[ind, ]) +
  geom_point(aes(sd_af, sd_ss)) +
  geom_abline(color = "red", size = 2)
# ggsave("tmp-figures/qcplot1.png", width = 7, height = 6)

ggplot(df[ind, ]) +
  geom_point(aes(sd_af, sd_ss, color = INFO)) +
  scale_color_viridis_c()
# ggsave("tmp-figures/qcplot2.png", width = 7, height = 6)

ggplot(df[ind, ]) +
  geom_point(aes(sd_af, sd_ss / sqrt(INFO), color = INFO)) +
  scale_color_viridis_c()
# ggsave("tmp-figures/qcplot3.png", width = 7, height = 6)

ggplot(df[ind, ]) +
  geom_point(aes(sd_af, sd_ss / sqrt(INFO),
                 color = ifelse(chr %in% c(6, 8, 11), chr, "Other"))) +
  labs(color = "Chr")
# ggsave("tmp-figures/qcplot4.png", width = 7, height = 6)

library(dplyr)
df %>%
  filter((sd_ss / sqrt(INFO) / sd_af) < 2.5) %>%
  select(chr, pos) %>%
  qplot(y = pos, color = as.factor(chr), data = .) +
  scale_y_continuous(breaks = 0:20 * 10e6)
