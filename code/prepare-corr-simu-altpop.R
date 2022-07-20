library(dplyr)
library(ggplot2)
library(bigsnpr)
library(bigreadr)
NCORES <- nb_cores()

snp_id <- with(snp_attach("data/ukbb4simu.rds")$map,
               paste(chromosome, physical.pos, allele1, allele2, sep = "_"))

code_country <- filter(fread2("UKBB/coding89.tsv"), selectable == "Y")
df0 <- fread2(
  "UKBB/ukb41181.csv",
  select    = c("eid", "20115-0.0", paste0("22009-0.", 1:16)),
  col.names = c("eid", "country",   paste0("PC", 1:16))
) %>%
  mutate(country = factor(country, levels = code_country$coding,
                          labels = code_country$meaning))

# still in data
sample <- fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
ind.indiv <- match(df0$eid, sample$ID_2)
id_to_rm <- scan("../../UKBiobank/ukb_withdrawn_individuals_8Feb2021.id")
ind.indiv[df0$eid %in% id_to_rm] <- NA
sub <- which(!is.na(ind.indiv))

# filter for relatedness
rel <- fread2("UKBB/ukb58024_rel_s488264.dat")
rel2 <- dplyr::filter(rel, Kinship > 2^-3.5)
ind.indiv[df0$eid %in% rel2$ID2] <- NA
sub <- which(!is.na(ind.indiv))

# Closer to South Europe (than UK)
PC_UKBB <- select(df0, PC1:PC16)[sub, ]
all_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)
center_italy <- unlist(all_centers[all_centers$Ancestry == "Italy", -1])
all_sq_dist <- rowSums(sweep(PC_UKBB, 2, center_italy, '-')^2)
closest_ind <- head(order(all_sq_dist), 10e3)
hist(log(all_sq_dist), "FD", xlim = c(6, 9)); abline(v = log(all_sq_dist[tail(closest_ind, 1)]), col = "red")
sub2 <- sub[closest_ind]
length(sub2) # 10,000
eid_UK <- df0$eid[sub2]

print_table <- function(x, min_count = 10) {
  tab <- sort(table(x, exclude = NULL), decreasing = TRUE)
  tab <- tab[tab >= min_count]
  cat(paste0(names(tab), ": ", tab), sep = " || ")
}
print_table(df0$country[sub2])
# NA: 6068 || Italy: 757 || France: 569 || Spain: 342 || Portugal: 279 || Germany: 187 ||
# USA: 146 || Cyprus: 143 || Malta: 132 || Switzerland: 127 || Greece: 117 || Brazil: 98 ||
# Egypt: 88 || South Africa: 69 || Bulgaria: 66 || Australia: 54 || Canada: 51 || Romania: 49 ||
# Austria: 46 || Serbia/Montenegro: 43 || Turkey: 37 || Argentina: 35 || Gibraltar: 34 ||
# Zimbabwe: 23 || Belgium: 23 || Croatia: 21 || Caribbean: 18 || Macedonia: 18 || Algeria: 17 ||
# Morocco: 16 || Bosnia and Herzegovina: 16 || United Kingdom: 15 || India: 14 || Albania: 12 ||
# Libya: 10 || Tunisia: 10 || Iraq: 10 || Hungary: 10

# Read the data
snp_readBGEN(
  bgenfiles = "UKBB/bgen/ukb_imp_chr22_v3.bgen",
  list_snp_id = list(snp_id),
  backingfile = "data/ukbb4simu_altpop",
  ind_row = match(df0$eid[sub2], sample$ID_2),
  read_as = "random",
  ncores = NCORES
)


library(bigsnpr)
ukb0 <- snp_attach("data/ukbb4simu.rds")
ukb <- snp_attach("data/ukbb4simu_altpop.rds")

snp_fst(list(
  data.frame(af = ukb0$map$freq, N = nrow(ukb0$genotypes)),
  data.frame(af =  ukb$map$freq, N = nrow( ukb$genotypes))
), overall = TRUE) # 0.001083306

qplot(ukb0$map$freq, ukb$map$freq, alpha = I(0.2)) +
  theme_bw(14) +
  coord_equal() +
  geom_abline(color = "red", linetype = 2) +
  labs(x = "Allele frequencies in the N.W. European LD reference",
       y = "Allele frequencies in the S. European LD reference")
# ggsave("figures/af_alt_pop.png", width = 7, height = 7)


#### Normal correlation matrix used in LDpred2 ####

G <- ukb$genotypes
POS2 <- snp_asGeneticPos(as.integer(ukb$map$chromosome),
                         ukb$map$physical.pos, dir = "tmp-data")

corr0 <- runonce::save_run(
  snp_cor(G, infos.pos = POS2, size = 3 / 1000, ncores = NCORES),
  file = "tmp-data/corr0_simu_altpop.rds")
as_SFBM(corr0, "tmp-data/corr_simu_altpop", compact = TRUE)$save()

corr0_NW <- readRDS("tmp-data/corr0_simu_val.rds")
ind <- sample(Matrix::which(corr0_NW != 0 & corr0 != 0), 100e3)
qplot(corr0_NW[ind], corr0[ind]) +
  theme_bw(14) +
  coord_equal() +
  geom_abline(color = "red", linetype = 2) +
  labs(x = "Pairwise correlations in the N.W. European LD reference",
       y = "Pairwise correlations in the S. European LD reference")
# ggsave("figures/corr_alt_pop.png", width = 7, height = 7)


#### Further restrict to nearly-independent blocks ####

splits <- snp_ldsplit(corr0, thr_r2 = 0.02, min_size = 50, max_size = 3:10 * 1000,
                      max_K = 50, max_r2 = 0.2, max_cost = 100)

splits$cost2 <- sapply(splits$all_size, function(sizes) sum(sizes^2))

library(dplyr)
best_split <- splits %>%
  arrange(cost2 * sqrt(10 + cost)) %>%
  print() %>%
  slice(1) %>%
  print()
# max_size n_block  cost perc_kept all_last   all_size       cost2
#     6000      12  10.5     0.716 <int [12]> <dbl [12]> 156225290
(all_size <- best_split$all_size[[1]])
# 2460 1554 2460 4806 3407 1685 3692 4502 5877 4911 2961 1685

final_grp <- rep(seq_along(all_size), all_size)

corr0T <- as(corr0, "dgTMatrix")
corr0T@x <- ifelse(final_grp[corr0T@i + 1L] == final_grp[corr0T@j + 1L], corr0T@x, 0)

as_SFBM(Matrix::drop0(corr0T), "tmp-data/corr_simu_altpop_with_blocks", compact = TRUE)$save()


#### Prepare per-block LD matrices for PRS-CS ####

bigassertr::assert_dir("tmp-data/ldblk_ukbb_simu_altpop")
hdf5_file <- "tmp-data/ldblk_ukbb_simu_altpop/ldblk_ukbb_chr22.hdf5"

ind_block <- split(seq_along(final_grp), final_grp)

for (ic in seq_along(ind_block)) {

  print(ic)
  ind <- ind_block[[ic]]
  ld <- big_cor(G, ind.col = ind)[]

  rhdf5::h5write(list(ldblk = ld, snplist = ukb$map$marker.ID[ind]),
                 file = hdf5_file, name = paste0("blk_", ic), level = 9)
}

# need also some .bim and some snpinfo (not hm3, but this is the filename expected)
ukb$map %>%
  mutate(allele1 = "A", allele2 = "C", genetic.dist = 0) %>%
  select(bigsnpr:::NAMES.MAP) %>%
  bigsnpr:::write.table2("tmp-data/ldblk_ukbb_simu_altpop/for_prscs.bim")

snpinfo <- transmute(ukb$map, CHR = as.integer(chromosome), SNP = marker.ID,
                     BP = physical.pos, A1 = "A", A2 = "C")
snpinfo$MAF <- snp_MAF(G, ncores = NCORES)
bigreadr::fwrite2(snpinfo, "tmp-data/ldblk_ukbb_simu_altpop/snpinfo_ukbb_hm3", sep = "\t")


#### Prepare shrunk LD matrix from GCTB ####

ukb$map$genetic.dist <- POS2
ukb$fam <- snp_fake(n = nrow(G), m = 1)$fam
bedfile <- "tmp-data/for_gctb_simu_altpop.bed"
snp_writeBed(ukb, bedfile)

gctb <- "tmp-data/gctb_2.03beta_Linux/gctb"

system(glue::glue(
  "{gctb} --bfile {sub_bed(bedfile)}",
  " --make-shrunk-ldm",
  " --out tmp-data/for_gctb_simu_altpop"
)) # Computational time: 1:27:14

system(glue::glue(
  "{gctb} --ldm tmp-data/for_gctb_simu_altpop.ldm.shrunk",
  " --make-sparse-ldm --chisq 0",
  " --out tmp-data/for_gctb_simu_altpop"
))  # Computational time: 0:6:40

# Convert their LD matrix to R

binfile <- "tmp-data/for_gctb_simu_altpop.ldm.sparse.bin"
info <- bigreadr::fread2(sub("\\.bin$", ".info", binfile))
I <- list()
X <- list()
con <- file(binfile, open = "rb")
for (n in info$WindSize) {
  cat(".")
  I[[length(I) + 1]] <- readBin(con, what = 1L, size = 4, n = n) + 1L
  X[[length(X) + 1]] <- readBin(con, what = 1,  size = 4, n = n)
}
close(con)

corr3 <- Matrix::sparseMatrix(i = unlist(I), j = rep(seq_along(I), lengths(I)), x = unlist(X))
corr3 <- Matrix::cov2cor(corr3)  # force diag to be exactly 1.0
as_SFBM(corr3, backingfile = "tmp-data/for_gctb_simu_altpop", compact = TRUE)$save()

# Add blocks to the shrunk LD matrix
corr3T <- as(corr3, "dgTMatrix")
corr3T@x <- ifelse(final_grp[corr3T@i + 1L] == final_grp[corr3T@j + 1L], corr3T@x, 0)

as_SFBM(Matrix::drop0(corr3T), "tmp-data/for_gctb_simu_altpop_with_blocks", compact = TRUE)$save()

file.size("tmp-data/for_gctb_simu_altpop.sbk") /
  file.size("tmp-data/for_gctb_simu_altpop_with_blocks.sbk") # 2.9
