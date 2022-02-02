library(bigsnpr)
ukb <- snp_attach("data/UKBB_large_val.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)
POS <- ukb$map$physical.pos
INFO <- ukb$map$info


#### Compute LD matrices ####

bigassertr::assert_dir("data/corr-large")

library(future.batchtools)
NCORES <- 30
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "200g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

future.apply::future_lapply(1:22, function(chr) {

  ind.chr <- which(CHR == chr)
  POS2 <- snp_asGeneticPos(CHR[ind.chr], POS[ind.chr], dir = "tmp-data")

  corr <- snp_cor(G, ind.col = ind.chr, infos.pos = POS2, size = 3 / 1000,
                  info = INFO[ind.chr], thr_r2 = 0.01, ncores = NCORES)

  saveRDS(corr, file = paste0("data/corr-large/chr", chr, ".rds"), version = 2)
})
# 10 min for chr 22

sum(sapply(list.files("data/corr-large", full.names = TRUE), file.size)) / 1024^3
# 41.7 GB
