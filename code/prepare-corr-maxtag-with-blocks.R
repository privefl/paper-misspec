library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = 2, mem = "150g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

all_final_grp <- furrr::future_map_dfr(1:22, function(chr) {

  corr0 <- readRDS(paste0("data/corr-large/maxtag_adj_chr", chr, ".rds"))
  dim(corr0)

  library(bigsnpr)
  all_splits <- runonce::save_run({
    SEQ <- round(seq_log(1000 + ncol(corr0) / 50,
                         5000 + ncol(corr0) / 10,
                         length.out = 10))

    splits <- purrr::map_dfr(SEQ, function(max_size) {
      res <- snp_ldsplit(corr0, thr_r2 = 0.02, min_size = 150,
                         max_size = max_size, max_K = 250)
      res$max_size <- max_size
      res
    })
  }, file = paste0("tmp-data/split_corr_maxtag_chr", chr, ".rds"))


  library(ggplot2)
  qplot(data = all_splits, perc_kept, cost, color = as.factor(max_size)) +
    theme_bw(12) +
    theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
    labs(y = "Sum of squared correlations outside blocks (log-scale)",
         x = "% of non-zero values kept", color = "Maximum block size") +
    ylim(0, min(quantile(all_splits$cost, 0.6), 200)) +
    geom_vline(xintercept = 0.6, linetype = 3)


  library(dplyr)
  final_split <- all_splits %>%
    filter(perc_kept < `if`(chr == 6, 0.7, 0.6)) %>%
    arrange(cost) %>%
    slice(1)

  final_grp <- final_split$block_num[[1]]

  corr0T <- as(corr0, "dgTMatrix")
  corr0T@x <- ifelse(final_grp[corr0T@i + 1L] == final_grp[corr0T@j + 1L], corr0T@x, 0)

  corr0T %>%
    Matrix::drop0() %>%
    as("symmetricMatrix") %>%
    saveRDS(paste0("data/corr-large/maxtag_adj_with_blocks_chr", chr, ".rds"))

  final_split
})

plot(all_final_grp$n_block)
plot(all_final_grp$cost)

sum(file.size(paste0("data/corr-large/maxtag_adj_with_blocks_chr", 1:22, ".rds"))) /
  sum(file.size(paste0("data/corr-large/maxtag_adj_chr", 1:22, ".rds")))
# 60.5%
