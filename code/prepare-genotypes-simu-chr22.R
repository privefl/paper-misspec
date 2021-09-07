# Get variant info
chr <- 22
info_chr <- bigreadr::fread2(paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt"))
info_chr_sub <- subset(info_chr, V6 > 0.01 & V8 > 0.4, -c(V1, V2, V7))
map <- cbind.data.frame(chr = chr, setNames(info_chr_sub, c("pos", "a1", "a2", "maf", "info")))
hist(map$info)
M <- 40e3

{
  all_info <- map$info
  kern_smooth <- ks::kde(all_info)
  set.seed(1)
  ind_snp <- sort(sample(length(all_info), M, replace = FALSE,
                         prob = 1 / predict(kern_smooth, x = all_info)))
  # pdf("figures/simu_hist_info.pdf", width = 7, height = 5)
  hist(all_info[ind_snp], main = NULL, xlab = "INFO")
  # dev.off()
  list_snp_id <- list(with(map[ind_snp, ], paste(chr, pos, a1, a2, sep = "_")))
}


df0 <- bigreadr::fread2(
  "UKBB/ukb41181.csv",
  select = c("eid", "22020-0.0", paste0("22009-0.", 1:16)),
  col.names = c("eid", "used_in_pca", paste0("PC", 1:16))
)

# sample still in data
sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
ind.indiv <- match(df0$eid, sample$ID_2)
id_to_rm <- scan("../../UKBiobank/ukb_withdrawn_individuals_8Feb2021.id")
ind.indiv[df0$eid %in% id_to_rm] <- NA
sub <- which(!is.na(ind.indiv) & df0$used_in_pca)

# Genetically homogeneous
dist <- bigutilsr::dist_ogk(as.matrix(df0[sub, -(1:2)]))
hist(log(dist), "FD")
sub2 <- sub[log(dist) < 5]
length(sub2) # 362,307


library(bigsnpr)
NCORES <- nb_cores()

# Read all data as random hard calls
bgen <- "UKBB/bgen/ukb_imp_chr22_v3.bgen"
system.time(
  snp_readBGEN(
    bgenfiles = bgen,
    list_snp_id = list_snp_id,
    backingfile = "data/ukbb4simu",
    ind_row = ind.indiv[sub2],
    ncores = NCORES,
    read_as = "random"
  )
) # 12 min

simu <- snp_attach("data/ukbb4simu.rds")
simu$map$info <- all_info[ind_snp]
simu$fam <- data.frame(eid = df0$eid[sub2], ind_bgen = ind.indiv[sub2])
snp_save(simu)
G <- simu$genotypes
G$file_size / 1024^3  # 13.5 GB
G[, 1]
G[, ncol(G)]

set.seed(1)
ind.val <- sort(sample(nrow(G), 10e3))
ind.gwas <- sort(sample(setdiff(rows_along(G), ind.val), 300e3))
ind.test <- setdiff(rows_along(G), c(ind.val, ind.gwas))
save(ind.gwas, ind.val, ind.test, file = "data/ukbb4simu_ind.RData")

# Read data as dosages
snp_readBGEN(
  bgenfiles = bgen,
  list_snp_id = list_snp_id,
  backingfile = "data/ukbb4simu_imp",
  ind_row = ind.indiv[sub2],
  ncores = NCORES
)
