library(dplyr)
library(ggplot2)
library(bigsnpr)
NCORES <- nb_cores()

snp_id <- with(snp_attach("data/ukbb4simu.rds")$map,
               paste(chromosome, physical.pos, allele1, allele2, sep = "_"))

df0 <- bigreadr::fread2(
  "UKBB/ukb41181.csv",
  select    = c("eid", paste0("22009-0.", 1:16)),
  col.names = c("eid", paste0("PC", 1:16))
)

# still in data
sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
ind.indiv <- match(df0$eid, sample$ID_2)
id_to_rm <- scan("../../UKBiobank/ukb_withdrawn_individuals_8Feb2021.id")
ind.indiv[df0$eid %in% id_to_rm] <- NA
sub <- which(!is.na(ind.indiv))

# filter for relatedness
rel <- bigreadr::fread2("UKBB/ukb58024_rel_s488264.dat")
rel2 <- dplyr::filter(rel, Kinship > 2^-3.5)
ind.indiv[df0$eid %in% rel2$ID2] <- NA
sub <- which(!is.na(ind.indiv))

# Closer to South Europe (than UK)
PC_UKBB <- select(df0, PC1:PC16)[sub, ]
all_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)
center_italy <- unlist(all_centers[4, -1])
all_sq_dist <- rowSums(sweep(PC_UKBB, 2, center_italy, '-')^2)
closest_ind <- head(order(all_sq_dist), 10e3)
hist(log(all_sq_dist), "FD", xlim = c(6, 9)); abline(v = log(all_sq_dist[tail(closest_ind, 1)]), col = "red")
sub2 <- sub[closest_ind]
length(sub2) # 10,000
eid_UK <- df0$eid[sub2]

# Read the data
snp_readBGEN(
  bgenfiles = "UKBB/bgen/ukb_imp_chr22_v3.bgen",
  list_snp_id = list(snp_id),
  backingfile = "data/ukbb4simu_altpop",
  ind_row = match(df0$eid[sub2], sample$ID_2),
  ncores = NCORES
)


library(bigsnpr)
ukb0 <- snp_attach("data/ukbb4simu.rds")
ukb <- snp_attach("data/ukbb4simu_altpop.rds")

snp_fst(list(
  data.frame(af = ukb0$map$freq, N = nrow(ukb0$genotypes)),
  data.frame(af =  ukb$map$freq, N = nrow( ukb$genotypes))
), overall = TRUE) # 0.001083306


#### Normal correlation matrix used in LDpred2 ####

G <- ukb$genotypes
POS2 <- snp_asGeneticPos(as.integer(ukb$map$chromosome),
                         ukb$map$physical.pos, dir = "tmp-data")

corr0 <- runonce::save_run(
  snp_cor(G, infos.pos = POS2, size = 3 / 1000, ncores = NCORES),
  file = "tmp-data/corr0_simu_altpop.rds")

corr <- runonce::save_run(as_SFBM(corr0, "tmp-data/corr_simu_altpop", compact = TRUE),
                          file = "tmp-data/corr_simu_altpop.rds")

corr_UK <- readRDS("tmp-data/corr0_simu_val.rds")
ind <- sample(Matrix::which(corr0 != 0 & corr_UK != 0), 100e3)
qplot(corr0[ind], corr_UK[ind]) + theme_bw() + geom_abline(color = "red")

#### Further restrict to nearly-independent blocks ####

splits <- purrr::map_dfr(3:8 * 1000, function(max_size) {
  res <- snp_ldsplit(corr0, thr_r2 = 0.02, min_size = 50,
                     max_size = max_size, max_K = 50)
  res$max_size <- max_size
  res
})

qplot(data = splits, perc_kept, cost, color = as.factor(max_size)) +
  theme_bw(12) +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
  labs(y = "Sum of squared correlations outside blocks (log-scale)",
       x = "% of non-zero values kept", color = "Maximum block size") +
  ylim(0, median(splits$cost)) +
  geom_vline(xintercept = 0.6, linetype = 3)

final_grp <- splits %>%
  filter(perc_kept < 0.6) %>%
  arrange(cost) %>%
  slice(1) %>%
  print() %>%
  pull(block_num) %>%
  unlist()

corr0T <- as(corr0, "dgTMatrix")
corr0T@x <- ifelse(final_grp[corr0T@i + 1L] == final_grp[corr0T@j + 1L], corr0T@x, 0)

corr <- runonce::save_run(
  as_SFBM(Matrix::drop0(corr0T), "tmp-data/corr_simu_altpop_with_blocks", compact = TRUE),
  file = "tmp-data/corr_simu_altpop_with_blocks.rds")
