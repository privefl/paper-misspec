library(bigsnpr)
ukb <- snp_attach("data/ukbb4simu.rds")

#### Normal correlation matrix used in LDpred2 ####

G <- ukb$genotypes
POS2 <- snp_asGeneticPos(as.integer(ukb$map$chromosome),
                         ukb$map$physical.pos, dir = "tmp-data")

load("data/ukbb4simu_ind.RData")

corr0 <- runonce::save_run(
  snp_cor(G, ind.row = ind.val, infos.pos = POS2, size = 3 / 1000, ncores = nb_cores()),
  file = "tmp-data/corr0_simu_val.rds")

corr <- runonce::save_run(as_SFBM(corr0, "tmp-data/corr_simu_val"),
                          file = "tmp-data/corr_simu_val.rds")


#### Further restrict to nearly-independent blocks ####

splits <- purrr::map_dfr(3:8 * 1000, function(max_size) {
  res <- snp_ldsplit(corr0, thr_r2 = 0.02, min_size = 50,
                     max_size = max_size, max_K = 50)
  res$max_size <- max_size
  res
})

corr0T <- as(corr0, "dgTMatrix")
corr0_i <- corr0T@i + 1L
corr0_j <- corr0T@j + 1L

library(future.apply)
plan("multisession", workers = nb_cores())
options(future.globals.maxSize = 2e9)
splits$cost_size <- future_sapply(splits$block_num, function(grp) {
  sum(grp[corr0_i] == grp[corr0_j])
})
plan("sequential")

library(ggplot2)
qplot(data = splits, cost_size / length(corr0T@x), cost, color = as.factor(max_size)) +
  theme_bw(12) +
  theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
  labs(y = "Sum of squared correlations outside blocks (log-scale)",
       x = "% of non-zero values kept", color = "Maximum block size") +
  ylim(0, median(splits$cost)) +
  geom_vline(xintercept = 0.6, linetype = 3)
# ggsave("figures/cost_split_chr22.pdf", width = 7, height = 6)

library(dplyr)
final_grp <- splits %>%
  filter(cost_size < (0.6 * length(corr0T@x))) %>%
  arrange(cost) %>%
  slice(1) %>%
  print() %>%
  pull(block_num) %>%
  unlist()

corr0T@x <- ifelse(final_grp[corr0_i] == final_grp[corr0_j], corr0T@x, 0)

corr <- runonce::save_run(
  as_SFBM(Matrix::drop0(corr0T), "tmp-data/corr_simu_val_with_blocks"),
  file = "tmp-data/corr_simu_val_with_blocks.rds")
