library(bigsnpr)
ukb <- snp_attach("data/UKBB_HM3_val.rds")
G <- ukb$genotypes
CHR <- ukb$map$chromosome
POS2 <- ukb$map$genetic.dist

(NCORES <- nb_cores())

bigassertr::assert_dir("data/corr")

for (chr in 1:22) {

  ind.chr <- which(CHR == chr)

  corr <- runonce::save_run(
    snp_cor(G, ind.col = ind.chr, infos.pos = POS2[ind.chr],
            size = 3 / 1000, ncores = NCORES),
    file = paste0("data/corr/chr", chr, ".rds")
  )
}
