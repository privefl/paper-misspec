library(bigsnpr)
library(furrr)
plan("multisession", workers = 6)

bigassertr::assert_dir("ldref")

all_final_grp <- future_map_dfr(1:22, function(chr) {

  corr0 <- readRDS(paste0("../paper-ldpred2/ld-ref/LD_chr", chr, ".rds"))
  dim(corr0)

  # find nearly independent LD blocks
  m <- ncol(corr0)
  (SEQ <- round(seq_log(m / 30, m / 5, length.out = 20)))
  splits <- snp_ldsplit(corr0, thr_r2 = 0.02, min_size = 50, max_size = SEQ, max_r2 = 0.1)

  library(dplyr)
  best_split <- splits %>%
    arrange(cost2 * (5 + cost)) %>%  # removed the sqrt() because costs were too high
    print() %>%
    slice(1) %>%
    print()

  (all_size <- best_split$all_size[[1]])
  best_grp <- rep(seq_along(all_size), all_size)

  corr0T <- as(corr0, "dgTMatrix")
  corr0T@x <- ifelse(best_grp[corr0T@i + 1L] == best_grp[corr0T@j + 1L], corr0T@x, 0)

  corr0T %>%
    Matrix::drop0() %>%
    as("symmetricMatrix") %>%
    saveRDS(file = paste0("ldref/LD_with_blocks_chr", chr, ".rds"), version = 2)

  best_split
}, .options = furrr_options(scheduling = FALSE))

plot(all_final_grp$n_block)
plot(all_final_grp$cost)
plot(all_final_grp$cost2)

sum(file.size(paste0("ldref/LD_with_blocks_chr", 1:22, ".rds"))) /
  sum(file.size(paste0("../paper-ldpred2/ld-ref/LD_chr", 1:22, ".rds")))
# 84%

# Add positions for different genome builds
liftOver <- runonce::download_file(
  "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver", "tmp-data")
map <- readRDS(url("https://ndownloader.figshare.com/files/25503788"))
sum(is.na(map))  # before: 709 (positions not matched)
map$pos_hg17 <- snp_modifyBuild(map, liftOver, from = "hg19", to = "hg17")$pos
map$pos_hg18 <- snp_modifyBuild(map, liftOver, from = "hg19", to = "hg18")$pos
map$pos_hg38 <- snp_modifyBuild(map, liftOver, from = "hg19", to = "hg38")$pos
sum(is.na(map))  # after: 909 (positions not matched with extra QC)

block_sizes <- unlist(all_final_grp$all_size)
map$group_id <- rep(seq_along(block_sizes), block_sizes)
saveRDS(map, "ldref/map.rds", version = 2)
system("zip -r -0 ldref_with_blocks.zip ldref/*")
