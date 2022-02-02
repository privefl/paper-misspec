library(bigsnpr)
ukb <- snp_attach("data/UKBB_large_val.rds")
G <- ukb$genotype
CHR <- ukb$map$chromosome <- as.integer(ukb$map$chromosome)
POS <- ukb$map$physical.pos

map_hm3 <- snp_attach("data/UKBB_HM3_val.rds")$map
in_hm3 <- vctrs::vec_in(ukb$map[1:6], map_hm3[1:6])
sum(in_hm3)  # 1,041,651

ind_keep <- runonce::save_run({

  exclude <- with(ukb$map, info < 0.4 | freq < 0.01 | freq > 0.99 |
                    paste(allele1, allele2) %in% c("A T", "T A", "C G", "G C"))

  snp_clumping(G, thr.r2 = 0.9, S = ukb$map$mean_info + in_hm3,
               infos.chr = CHR, infos.pos = POS, size = 100,
               exclude = which(exclude), ncores = nb_cores())
}, file = "data/ind_keep_large_clump.rds")
# 9 min -- 2.5 M variants with thr_r2 = 0.9

sum(in_hm3[ind_keep]) # 554,655
table(CHR_keep <- CHR[ind_keep])
#      1      2      3      4      5      6      7      8      9     10     11
# 190939 204351 170390 164476 151197 150797 146197 132861 112483 123023 119075
#     12     13     14     15     16     17     18     19     20     21     22
# 115771  88334  78674  76638  86605  78339  74165  65514  60271  35898  39480

POS2_keep <- snp_asGeneticPos(CHR_keep, POS[ind_keep], dir = "tmp-data",
                              ncores = nb_cores())
INFO_keep <- ukb$map$info[ind_keep]

bigassertr::assert_dir("data/corr-large")

library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_walk(1:22, function(chr) {

  ind.chr <- which(CHR_keep == chr)

  runonce::save_run(
    snp_cor(G, ind.col = ind_keep[ind.chr], infos.pos = POS2_keep[ind.chr],
            size = 3 / 1000, info = INFO_keep[ind.chr], ncores = NCORES),
    file = paste0("data/corr-large/clump_adj_chr", chr, ".rds")
  )
})
