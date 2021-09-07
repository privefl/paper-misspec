library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = 2, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_walk(1:22, function(chr) {

  runonce::save_run({

    corr2 <- as(readRDS(paste0("data/corr-large/chr", chr, ".rds")), "dgCMatrix")

    library(bigsnpr)
    ukbb <- snp_attach("data/UKBB_large_val.rds")
    map_large <- dplyr::mutate(ukbb$map, chromosome = as.integer(chromosome))
    ind.chr <- which(map_large$chromosome == chr)
    map_large.chr <- map_large[ind.chr, ]

    map_hm3 <- snp_attach("data/UKBB_HM3_val.rds")$map
    keep <- vctrs::vec_in(map_large.chr[1:6], map_hm3[1:6])

    MAF <- with(map_large.chr, pmin(freq, 1 - freq))
    exclude <- MAF < 0.01 | map_large.chr$info < 0.3 |
      with(map_large.chr, paste(allele1, allele2) %in% c("A T", "T A", "C G", "G C"))

    Rcpp::sourceCpp('code/greedy-maxtag.cpp')
    set_max_tag(corr2@p, corr2@i, corr2@x, min_add = 0.2,
                select = keep, exclude = exclude, remove_diag = FALSE,
                sqrt_info = sqrt(map_large.chr$mean_info))

  }, file = paste0("tmp-data/tag_set_chr", chr, ".rds"))
})
# 30 min max for all chromosomes but #6


## Prepare validation set
library(bigsnpr)
ukbb <- snp_attach("data/UKBB_large_val.rds")

res <- unlist(lapply(1:22, function(chr) {
  readRDS(paste0("tmp-data/tag_set_chr", chr, ".rds"))[[1]]
}))
stopifnot(length(res) == nrow(ukbb$map))
(ind <- which(!is.na(res)))  # 2.03M
ukbb$fam <- snp_attach("data/UKBB_HM3_val.rds")$fam
snp_subset(ukbb, ind.col = ind, backingfile = "data/UKBB_maxtag_val")


## Prepare test set
list_snp_id <- with(ukbb$map[ind, ],
                    split(paste(chromosome, physical.pos, allele1, allele2, sep = "_"),
                          as.integer(chromosome)))
# saveRDS(list_snp_id, "data/UKBB_maxtag_snp_id.rds")
