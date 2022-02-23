library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

map <- transmute(snp_attach("data/UKBB_HM3_val.rds")$map,
                 chr = as.integer(chromosome), pos = physical.pos, af_val = freq,
                 a0 = allele1, a1 = allele2)  # reversed somehow..


#### Type 2 diabetes (T2D) ####

# Download from http://diagram-consortium.org/downloads.html
# DIAGRAM 1000G GWAS meta-analysis Stage 1 Summary statistics
# Published in Scott et al (2017)
# unzip("tmp-data/METAANALYSIS_DIAGRAM_SE1.zip", exdir = "tmp-data")
sumstats <- fread2("tmp-data/METAANALYSIS_DIAGRAM_SE1.txt")
sumstats <- tidyr::separate(sumstats, "Chr:Position", c("chr", "pos"), convert = TRUE)
names(sumstats) <- c("chr", "pos", "a1", "a0", "beta", "beta_se", "p", "N")
head(sumstats$beta)     # -0.0130 -0.0710  0.0140  0.1500 -0.0085 -0.0060
head(sumstats$beta_se)  # 0.026 0.170 0.012 0.130 0.013 0.016

table_text <- "cases/408/347\ncontrols/3,098/3,911\ncases/97/35\ncontrols/297/158\ncases/181/74\ncontrols/801/846\ncases/3,982/3,357\ncontrols/35,946/47,103\ncases/413/266\ncontrols/281/416\ncases/530/493\ncontrols/541/533\ncases/41/39\ncontrols/866/902\ncases/161/228\ncontrols/2,754/3,259\ncases/2,229/2,395\ncontrols/1,673/2,995\ncases/386/287\ncontrols/3,441/4,219\ncases/653/508\ncontrols/574/600\ncases/1,773/1,525\ncontrols/1,376/1,267\ncases/-/1,124\ncontrols/-/1,298\ncases/201/146\ncontrols/764/786\ncases/-/1,467\ncontrols/-/1,754\ncases/65/46\ncontrols/409/429\ncases/251/403\ncontrols/2,111/3,108\ncases/166/0\ncontrols/953/0\ncases/1,118/806\ncontrols/1,446/1,492"
sizes <- bigreadr::fread2(text = gsub(",", "", table_text, fixed = TRUE),
                          header = FALSE, na = "-", sep = "/")

sum_sizes <- rowSums(sizes[2:3], na.rm = TRUE)  # makes + females
sum(Nca <- sum_sizes[seq(1, length(sum_sizes), by = 2)])  #  26,201
sum(Nco <- sum_sizes[seq(2, length(sum_sizes), by = 2)])  # 132,407
(Neff <- sum(4 / (1 / Nca + 1 / Nco)))  # 70872.33
4 / (1 / sum(Nca) + 1 / sum(Nco))  # 87491.07
max(sumstats$N)  # 158,186
26201 + 132407   # 158,608


info_snp <- bigsnpr::snp_match(sumstats, map)

info_snp2 <- info_snp %>%
  filter(N > (0.7 * max(N))) %>%
  mutate(
    n_eff = N * Neff / max(N),
    sd_af = sqrt(2 * af_val * (1 - af_val)),
    sd_ss = 2 / sqrt(n_eff * beta_se^2 + beta^2))

qplot(sd_af, sd_ss, alpha = I(0.1),
      data = slice_sample(info_snp2, n = 50e3)) +
  theme_bw(15) +
  coord_equal() +
  # scale_color_viridis_c(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies",
       y = "Standard deviations derived from the summary statistics")
