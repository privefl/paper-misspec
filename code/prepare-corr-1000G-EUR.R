library(dplyr)
library(bigreadr)
library(bigsnpr)

obj.1000G <- snp_attach("data/1000G_HM3.rds")
G <- obj.1000G$genotypes
ind_EUR <- which(obj.1000G$fam$`Super Population` == "EUR")


#### Prepare correlations ####

CHR <- obj.1000G$map$chromosome
POS2 <- snp_attach("data/UKBB_HM3_val.rds")$map$genetic.dist

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = 15, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("data/corr_1000G_EUR")

all_final_grp <- furrr::future_map_dfr(1:22, function(chr) {

  ind.chr <- which(CHR == chr)

  library(bigsnpr)
  corr0 <- runonce::save_run({
    corr0 <- snp_cor(G, ind.row = ind_EUR, ind.col = ind.chr, infos.pos = POS2[ind.chr],
                     size = 3 / 1000, ncores = nb_cores())
    mean(is.na(corr0@x))  # 0.001 for chr 22
    corr0@x <- ifelse(is.na(corr0@x) | is.infinite(corr0@x), 0, corr0@x)
    corr0
  }, file = paste0("data/corr_1000G_EUR/chr", chr, ".rds"))

  library(furrr)
  plan("multisession", workers = 5)
  options(future.globals.maxSize = 10e9)

  all_splits <- runonce::save_run({
    SEQ <- round(seq_log(1000 + ncol(corr0) / 50,
                         5000 + ncol(corr0) / 10,
                         length.out = 10))

    splits <- future_map_dfr(SEQ, function(max_size) {
      res <- snp_ldsplit(corr0, thr_r2 = 0.05, min_size = 100,
                         max_size = max_size, max_K = 200)
      res$max_size <- max_size
      res
    })
  }, file = paste0("tmp-data/split_corr_1000G_EUR_chr", chr, ".rds"))


  library(ggplot2)
  qplot(data = all_splits, perc_kept, cost, color = as.factor(max_size)) +
    theme_bw(12) +
    theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
    labs(y = "Sum of squared correlations outside blocks (log-scale)",
         x = "% of non-zero values kept", color = "Maximum block size") +
    ylim(0, min(quantile(all_splits$cost, 0.6), 2000)) +
    geom_vline(xintercept = 0.6, linetype = 3) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2))


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
    saveRDS(paste0("data/corr_1000G_EUR/adj_with_blocks_chr", chr, ".rds"))

  final_split
})

plot(all_final_grp$n_block)
plot(all_final_grp$cost)
