library(dplyr)
library(bigreadr)
library(bigsnpr)

obj.1000G <- snp_attach("data/1000G_HM3.rds")
G <- obj.1000G$genotypes
CHR <- obj.1000G$map$chromosome
ind_EUR <- which(obj.1000G$fam$`Super Population` == "EUR")


#### Prepare correlations ####

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = 15, mem = "60g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("data/corr_1000G_EUR")
bigassertr::assert_dir("data/corr_1000G_EUR/1kg_for_prscs")

all_final_grp <- furrr::future_map_dfr(1:22, function(chr) {

  ind.chr <- which(CHR == chr)
  MAF <- snp_MAF(G, ind.row = ind_EUR, ind.col = ind.chr, ncores = nb_cores())
  keep <- (MAF > 0.02)
  ind.chr2 <- ind.chr[keep]

  library(bigsnpr)
  corr0 <- runonce::save_run({
    POS2 <- snp_attach("data/UKBB_HM3_val.rds")$map$genetic.dist
    snp_cor(G, ind.row = ind_EUR, ind.col = ind.chr2, infos.pos = POS2[ind.chr2],
            size = 3 / 1000, ncores = nb_cores())
  }, file = paste0("data/corr_1000G_EUR/chr", chr, ".rds"))


  # find nearly independent LD blocks
  m <- length(ind.chr2)
  (SEQ <- round(seq_log(m / 30, m / 5, length.out = 20)))
  splits <- snp_ldsplit(corr0, thr_r2 = 0.05, min_size = 50, max_size = SEQ, max_r2 = 0.15)
  splits$cost2 <- sapply(splits$all_size, function(sizes) sum(sizes^2))

  best_split <- splits %>%
    arrange(cost2 * sqrt(5 + cost)) %>%
    print() %>%
    slice(1) %>%
    print()

  (all_size <- best_split$all_size[[1]])
  best_grp <- rep(seq_along(all_size), all_size)

  runonce::save_run({
    corr0T <- as(corr0, "dgTMatrix")
    corr0T@x <- ifelse(best_grp[corr0T@i + 1L] == best_grp[corr0T@j + 1L], corr0T@x, 0)
    as(Matrix::drop0(corr0T), "symmetricMatrix")
  }, file = paste0("data/corr_1000G_EUR/adj_with_blocks_chr", chr, ".rds"))


  # LD for PRS-CS (based on the same data and LD blocks)
  hdf5_file <- paste0("data/corr_1000G_EUR/1kg_for_prscs/ldblk_1kg_chr", chr, ".hdf5")
  runonce::skip_run_if({
    ind_block <- split(ind.chr2, best_grp)
    bigparallelr::set_blas_ncores(nb_cores())
    for (ic in seq_along(ind_block)) {
      ind <- ind_block[[ic]]
      ld <- big_cor(G, ind.row = ind_EUR, ind.col = ind)[]
      rhdf5::h5write(list(ldblk = ld, snplist = paste0("snp", ind)),
                     file = hdf5_file, name = paste0("blk_", ic), level = 9)
    }
  }, files = hdf5_file)


  # return
  tibble(best_split, ind = list(ind.chr2))
})


# verif
plot(all_final_grp$n_block)
plot(all_final_grp$cost)
plot(all_final_grp$cost2)
# saveRDS(all_final_grp, "data/corr_1000G_EUR/all_final_grp.rds")

sum(file.size(paste0("data/corr_1000G_EUR/adj_with_blocks_chr", 1:22, ".rds"))) /
  sum(file.size(paste0("data/corr_1000G_EUR/chr", 1:22, ".rds")))
# 63.0%


# need also some .bim and some snpinfo for PRS-CS
ind_keep <- unlist(all_final_grp$ind)
snpinfo <- transmute(obj.1000G$map, chromosome, marker.ID = paste0("snp", row_number()),
                     genetic.dist = 0, physical.pos, allele1 = "A", allele2 = "C")[ind_keep, ]
bigsnpr:::write.table2(snpinfo, "data/corr_1000G_EUR/1kg_for_prscs/for_prscs.bim")

snpinfo2 <- transmute(snpinfo, CHR = chromosome, SNP = marker.ID,
                      BP = physical.pos, A1 = allele1, A2 = allele2)
snpinfo2$MAF <- snp_MAF(G, ind.row = ind_EUR, ind.col = ind_keep, ncores = nb_cores())
bigreadr::fwrite2(snpinfo2, "data/corr_1000G_EUR/1kg_for_prscs/snpinfo_1kg_hm3", sep = "\t")
