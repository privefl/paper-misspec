library(bigsnpr)
ukb_test <- snp_attach("data/UKBB_HM3_test.rds")
G_test <- ukb_test$genotypes
ind.test <- ukb_test$fam$id_csv

library(dplyr)
grid <- tidyr::expand_grid(
  qc = c("noqc", "qc1", "qc2"),
  name_corr = c("data/corr/chr", "data/corr/adj_with_blocks_chr")
) %>%
  print()

library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("results/cad")

grid$res <- furrr::future_pmap(grid[1:2], function(qc, name_corr) {

  gwas_file <- paste0("data/sumstats/cad_", qc, ".rds")
  sumstats <- readRDS(gwas_file)
  num_id <- sumstats[["_NUM_ID_"]]

  basename <- file.path("results/cad", sub("\\.rds$", "", basename(gwas_file)))
  if (grepl("blocks", name_corr, fixed = TRUE))
    basename <- paste0(basename, "_adj_with_blocks")

  library(bigsnpr)
  CHR <- snp_attach("data/UKBB_HM3_val.rds")$map$chromosome

  tmp <- tempfile(tmpdir = "tmp-data/sfbm/")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

  for (chr in 1:22) {

    print(chr)

    ## indices in 'sumstats'
    ind.chr <- which(sumstats$chr == chr)
    ## indices in 'G'
    ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
    ## indices in 'corr'
    ind.chr3 <- match(ind.chr2, which(CHR == chr))

    corr0 <- readRDS(paste0(name_corr, chr, ".rds"))
    ld_chr <- Matrix::colSums(corr0^2)
    df_beta_chr <- sumstats[ind.chr, c("beta", "beta_se", "n_eff")]

    if (chr == 1) {
      df_beta <- df_beta_chr
      ld <- ld_chr[ind.chr3]
      corr <- as_SFBM(corr0[ind.chr3, ind.chr3], tmp, compact = TRUE)
    } else {
      df_beta <- rbind(df_beta, df_beta_chr)
      ld <- c(ld, ld_chr[ind.chr3])
      corr$add_columns(corr0[ind.chr3, ind.chr3], nrow(corr))
    }
  }


  # Run LD score regression
  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL)))
  h2_est <- ldsc[["h2"]]


  purrr::map_dfr(0:10 / 10, function(shrinkage) {

    prs_effects <- rep(0, ncol(G_test))
    prs_effects[num_id] <- runonce::save_run({
      multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est, ncores = NCORES,
                                     vec_p_init = seq_log(1e-4, 0.5, 30),
                                     allow_jump_sign = FALSE, shrink_corr = shrinkage)
      range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
      keep <- (range > (0.9 * quantile(range, 0.9)))
      rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    }, file = paste0(basename, "_", shrinkage, ".rds"))

    tibble::tibble(shrinkage = shrinkage, effects = list(prs_effects))
  })
})
# saveRDS(grid, "results/cad/all_res.rds")

grid2 <- tidyr::unnest(grid, "res")

bigparallelr::set_blas_ncores(NCORES)
all_pred <- big_prodMat(G_test, do.call("cbind", grid2$effects))
grid2$pred <- as.list(as.data.frame(all_pred))

df0 <- bigreadr::fread2(
  "UKBB/ukb41181.csv",
  select = c("eid", "34-0.0", "52-0.0", "22001-0.0", "21022-0.0", "189-0.0",
             paste0("22009-0.", 1:16)),
  col.names = c("eid", "year", "month", "sex", "age", "deprivation_index",
                paste0("PC", 1:16))
) %>%
  mutate(date = (year - 1900) + (month - 0.5) / 12, year = NULL, month = NULL)

covar.test <- as.matrix(df0[ind.test, -1])

y <- readRDS("data/all_phecodes.rds")[["411.4"]]
grid2$pcor <- purrr::map(grid2$pred, function(pred) {
  if (all(is.na(pred))) return(rep(NA, 3))
  pcor(pred, y[ind.test], covar.test)
})
# saveRDS(grid2, "results/cad/all_res2.rds")

all_res <- grid2 %>%
  mutate(use_blocks = grepl("with_blocks", name_corr)) %>%
  select(-effects, -pred, -name_corr) %>%
  tidyr::unnest_wider("pcor", names_sep = "_") %>%
  mutate(across(starts_with("pcor_"), function(x) sign(x) * x^2)) %>%
  print()


library(ggplot2)
ggplot(filter(all_res, !use_blocks), aes(as.factor(shrinkage), pcor_1, fill = qc)) +
  facet_wrap(~ I("cad"), ncol = 1) +
  bigstatsr::theme_bigstatsr(0.9) +
  scale_fill_manual(values =  c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = pcor_2, ymax = pcor_3),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  labs(x = "Shrinkage", y = "Partial phenotypic variance explained by PGS",
       fill = "Which QC?") +
  theme(legend.position = "top",
        legend.text = element_text(margin = margin(r = 10, unit = "pt")),
        legend.title = element_text(margin = margin(r = 10, unit = "pt"))) +
  geom_col(data = filter(all_res, use_blocks), position = position_dodge(),
           color = "red", alpha = 0, show.legend = FALSE)

# ggsave("figures/res-cad-shrinkage.pdf", width = 10, height = 6.5)
