library(bigsnpr)
ukb <- snp_attach("data/ukbb4simu.rds")

#### Normal correlation matrix used in LDpred2 ####

G <- ukb$genotypes
POS2 <- snp_asGeneticPos(as.integer(ukb$map$chromosome),
                         ukb$map$physical.pos, dir = "tmp-data")

load("data/ukbb4simu_ind.RData")

corr0 <- runonce::save_run(
  snp_cor(G, ind.row = ind.val, infos.pos = POS2, size = 3 / 1000, ncores = nb_cores()),
  file = "tmp-data/corr0_simu_val.rds")
print(range(corr0@x), digits = 22)

corr <- as_SFBM(corr0, "tmp-data/corr_simu_val", compact = TRUE)$save()


#### Further restrict to nearly-independent blocks ####

# /!\ need v1.10.1
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
#     6000      12  31.1     0.707 <int [12]> <dbl [12]> 154489062

library(ggplot2)
plot_grid(
  qplot(data = splits, perc_kept, cost, color = as.factor(max_size)) +
    geom_point(data = best_split, size = 2, color = "black") +
    theme_bw(12) +
    theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
    scale_x_continuous(limits = c(0, NA), breaks = seq(0, 1, by = 0.1)) +
    labs(y = "Sum of squared correlations outside blocks",
         x = "% of non-zero values kept", color = "Maximum block size"),
  qplot(data = splits, cost2, cost, color = as.factor(max_size)) +
    geom_point(data = best_split, size = 2, color = "black") +
    theme_bw(12) +
    theme(legend.position = "none") +
    xlim(0, NA) +
    labs(y = "Sum of squared correlations outside blocks",
         x = "Sum of squared blocks", color = "Maximum block size"),
  scale = 0.95, ncol = 1, rel_heights = c(1.18, 1)
)
# ggsave("figures/cost_split_chr22.pdf", width = 7, height = 8.5)


(all_size <- best_split$all_size[[1]])
# 2460 1554 2920 4346 3407 1685 3692 4501 5876 4913 2961 1685

final_grp <- rep(seq_along(all_size), all_size)

corr0T <- as(corr0, "dgTMatrix")
corr0T@x <- ifelse(final_grp[corr0T@i + 1L] == final_grp[corr0T@j + 1L], corr0T@x, 0)

corr2 <- as_SFBM(Matrix::drop0(corr0T), "tmp-data/corr_simu_val_with_blocks",
                 compact = TRUE)$save()


#### Prepare per-block LD matrices for PRS-CS ####

# plink <- download_plink("tmp-data")
# ukb$map$genetic.dist <- POS2
# ukb$fam <- snp_fake(n = nrow(G), m = 1)$fam

bigassertr::assert_dir("tmp-data/ldblk_ukbb_simu")
hdf5_file <- "tmp-data/ldblk_ukbb_simu/ldblk_ukbb_chr22.hdf5"

ind_block <- split(seq_along(final_grp), final_grp)

for (ic in seq_along(ind_block)) {

  print(ic)

  ind <- ind_block[[ic]]

  # tmp_bed <- snp_writeBed(ukb, tempfile(tmpdir = "tmp-data", fileext = ".bed"),
  #                         ind.row = ind.val, ind.col = ind)
  #
  # system(glue::glue(
  #   "{plink}",
  #   " --bfile {sub_bed(tmp_bed)}",
  #   " --keep-allele-order",
  #   " --r square",
  #   " --out tmp-data/corr_simu_prscs/part{ic}",
  #   " --memory 20000 --threads {nb_cores()}"
  # )) # 1 sec
  #
  # file.remove(paste0(sub_bed(tmp_bed), c(".bed", ".bim", ".fam")))
  #
  # ld <- as.matrix(bigreadr::fread2(paste0("tmp-data/corr_simu_prscs/part", ic, ".ld")))
  # colnames(ld) <- NULL

  # why use PLINK when you have bigsnpr? :')
  ld <- big_cor(G, ind.row = ind.val, ind.col = ind)[]

  # why use python when you have R? :')
  rhdf5::h5write(list(ldblk = ld, snplist = ukb$map$marker.ID[ind]),
                 file = hdf5_file, name = paste0("blk_", ic), level = 9)
}

# need also some .bim and some snpinfo (not hm3, but this is the filename expected)
ukb$map %>%
  mutate(allele1 = "A", allele2 = "C", genetic.dist = 0) %>%
  select(bigsnpr:::NAMES.MAP) %>%
  bigsnpr:::write.table2("tmp-data/ldblk_ukbb_simu/for_prscs.bim")

snpinfo <- transmute(ukb$map, CHR = as.integer(chromosome), SNP = marker.ID,
                     BP = physical.pos, A1 = "A", A2 = "C")
snpinfo$MAF <- snp_MAF(G, ind.row = ind.val, ncores = NCORES)
bigreadr::fwrite2(snpinfo, "tmp-data/ldblk_ukbb_simu/snpinfo_ukbb_hm3", sep = "\t")


#### Prepare shrunk LD matrix from GCTB ####

ukb$map$genetic.dist <- POS2
ukb$fam <- snp_fake(n = nrow(G), m = 1)$fam
bedfile <- "tmp-data/for_gctb_simu.bed"
snp_writeBed(ukb, bedfile, ind.row = ind.val)

unzip("tmp-data/gctb_2.03beta_Linux.zip", exdir = "tmp-data", overwrite = FALSE)
gctb <- "tmp-data/gctb_2.03beta_Linux/gctb"
bigsnpr:::make_executable(gctb)

system(glue::glue(
  "{gctb} --bfile {sub_bed(bedfile)}",
  " --make-shrunk-ldm",
  " --out tmp-data/for_gctb_simu"
)) # Computational time: 1:13:36

system(glue::glue(
  "{gctb} --ldm tmp-data/for_gctb_simu.ldm.shrunk",
  " --make-sparse-ldm --chisq 0",
  " --out tmp-data/for_gctb_simu"
))

# Convert their LD matrix to R

binfile <- "tmp-data/for_gctb_simu.ldm.sparse.bin"
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
as_SFBM(corr3, backingfile = "tmp-data/for_gctb_simu", compact = TRUE)$save()

# Add blocks to the shrunk LD matrix
corr3T <- as(corr3, "dgTMatrix")
corr3T@x <- ifelse(final_grp[corr3T@i + 1L] == final_grp[corr3T@j + 1L], corr3T@x, 0)

corr4 <- as_SFBM(Matrix::drop0(corr3T), "tmp-data/for_gctb_simu_with_blocks",
                 compact = TRUE)$save()

file.size("tmp-data/for_gctb_simu.sbk") /
  file.size("tmp-data/for_gctb_simu_with_blocks.sbk")  # 2.95
