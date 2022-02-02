library(bigsnpr)
ukb <- snp_attach("data/UKBB_maxtag_val.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)
table(CHR)
#      1      2      3      4      5      6      7      8      9     10     11
# 163457 168666 141715 133943 123370 128468 116418 110180  92967 103221 100611
#     12     13     14     15     16     17     18     19     20     21     22
#  95061  74998  64767  62391  66374  60870  62353  46982  51635  29977  30662

POS2 <- snp_asGeneticPos(CHR, ukb$map$physical.pos, dir = "tmp-data",
                         ncores = nb_cores())
INFO <- ukb$map$info

bigassertr::assert_dir("data/corr-large")

library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_walk(1:22, function(chr) {

  ind.chr <- which(CHR == chr)

  runonce::save_run(
    snp_cor(G, ind.col = ind.chr, infos.pos = POS2[ind.chr],
            size = 3 / 1000, info = INFO[ind.chr], ncores = NCORES),
    file = paste0("data/corr-large/maxtag_adj_chr", chr, ".rds")
  )
})
