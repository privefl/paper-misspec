library(bigsnpr)
library(dplyr)
ukb <- snp_attach("data/UKBB_HM3_val.rds")
G     <- ukb$genotypes
CHR   <- ukb$map$chromosome
POS2  <- ukb$map$genetic.dist

library(future.batchtools)
NCORES <- 14
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "50g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("data/corr")
bigassertr::assert_dir("data/corr/ukbb_for_prscs")

all_final_grp <- furrr::future_map_dfr(1:22, function(chr) {

  ind.chr <- which(CHR == chr)


  # LD for LDpred2 and lassosum2
  corr0 <- runonce::save_run(
    snp_cor(G, ind.col = ind.chr, infos.pos = POS2[ind.chr],
            size = 3 / 1000, ncores = NCORES),
    file = paste0("data/corr/chr", chr, ".rds")
  )


  # find nearly independent LD blocks
  m <- ncol(corr0)
  (SEQ <- round(seq_log(m / 30, m / 5, length.out = 20)))
  splits <- snp_ldsplit(corr0, thr_r2 = 0.02, min_size = 50, max_size = SEQ, max_r2 = 0.1)
  splits$cost2 <- sapply(splits$all_size, function(sizes) sum(sizes^2))

  best_split <- splits %>%
    arrange(cost2 * sqrt(10 + cost)) %>%
    print() %>%
    slice(1) %>%
    print()

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

  (all_size <- best_split$all_size[[1]])
  best_grp <- rep(seq_along(all_size), all_size)

  runonce::save_run({
    corr0T <- as(corr0, "dgTMatrix")
    corr0T@x <- ifelse(best_grp[corr0T@i + 1L] == best_grp[corr0T@j + 1L], corr0T@x, 0)
    as(Matrix::drop0(corr0T), "symmetricMatrix")
  }, file = paste0("data/corr/adj_with_blocks_chr", chr, ".rds"))


  # LD for PRS-CS (based on the same data and LD blocks)
  hdf5_file <- paste0("data/corr/ukbb_for_prscs/ldblk_ukbb_chr", chr, ".hdf5")
  runonce::skip_run_if({
    ind_block <- split(ind.chr, best_grp)
    bigparallelr::set_blas_ncores(NCORES)
    for (ic in seq_along(ind_block)) {
      ind <- ind_block[[ic]]
      ld <- big_cor(G, ind.col = ind)[]
      rhdf5::h5write(list(ldblk = ld, snplist = paste0("snp", ind)),
                     file = hdf5_file, name = paste0("blk_", ic), level = 9)
    }
  }, files = hdf5_file)


  # return
  best_split
})


# verif
plot(all_final_grp$n_block)
plot(all_final_grp$cost)
plot(all_final_grp$cost2)
# saveRDS(all_final_grp, "data/corr/all_final_grp.rds")

sum(file.size(paste0("data/corr/adj_with_blocks_chr", 1:22, ".rds"))) /
  sum(file.size(paste0("data/corr/chr", 1:22, ".rds")))
# 69.3%


# need also some .bim and some snpinfo for PRS-CS
snpinfo <- transmute(ukb$map, chromosome, marker.ID = paste0("snp", row_number()),
                     genetic.dist = 0, physical.pos, allele1 = "A", allele2 = "C")
bigsnpr:::write.table2(snpinfo, "data/corr/ukbb_for_prscs/for_prscs.bim")

snpinfo2 <- transmute(snpinfo, CHR = chromosome, SNP = marker.ID,
                      BP = physical.pos, A1 = allele1, A2 = allele2)
snpinfo2$MAF <- snp_MAF(G, ncores = NCORES)
bigreadr::fwrite2(snpinfo2, "data/corr/ukbb_for_prscs/snpinfo_ukbb_hm3", sep = "\t")
