# file.symlink("../../UKBB/", "UKBB")

map_hapmap3 <- readRDS(url("https://ndownloader.figshare.com/files/25503788"))

library(doParallel)
cl <- makeCluster(22)
parallel::clusterExport(cl, "map_hapmap3")
list_snp_id <- parLapply(cl, 1:22, function(chr) {
  mfi <- paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt")
  infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
  joined <- dplyr::inner_join(cbind(chr = chr, infos_chr), map_hapmap3[1:2],
                              by = c("chr", "V3" = "pos"))
  with(joined[!vctrs::vec_duplicate_detect(joined$V3), ],
       paste(chr, V3, V4, V5, sep = "_"))
})
stopCluster(cl)

sum(lengths(list_snp_id))  # 1,054,315


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
# saveRDS(ind.indiv[sub2], "tmp-data/eur_ind_bgen.rds")


# Validation set
set.seed(1); ind.val <- sort(sample(sub2, 10e3))

library(bigsnpr)
snp_readBGEN(
  bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
  list_snp_id = list_snp_id,
  backingfile = "data/UKBB_HM3_val",
  ind_row     = ind.indiv[ind.val],
  ncores      = nb_cores()
)

ukb <- snp_attach("data/UKBB_HM3_val.rds")
G <- ukb$genotypes
dim(G) # 10,000 x 1,054,315
G$file_size / 1024^3  # 10 GB
ukb$map$chromosome <- as.integer(ukb$map$chromosome)
ukb$map$genetic.dist <- snp_asGeneticPos(ukb$map$chromosome, ukb$map$physical.pos,
                                         dir = "../paper-ldpred2/tmp-data/",
                                         ncores = nb_cores())
str(ukb$map)

ukb$fam <- tibble::tibble(id_csv = ind.val)
snp_save(ukb)

# Test set
ind.test <- setdiff(sub2, ind.val)
snp_readBGEN(
  bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
  list_snp_id = list_snp_id,
  backingfile = "data/UKBB_HM3_test",
  ind_row     = ind.indiv[ind.test],
  ncores      = nb_cores()
)

ukb <- snp_attach("data/UKBB_HM3_test.rds")
G <- ukb$genotypes
dim(G) # 352,307 x 1,054,315
G$file_size / 1024^3  # 346 GB

ukb$fam <- tibble::tibble(id_csv = ind.test)
snp_save(ukb)
