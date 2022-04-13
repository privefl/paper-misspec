library(bigsnpr)
ukb <- snp_attach("data/UKBB_HM3_val.rds")
G <- ukb$genotypes
NCORES <- nb_cores()

library(dplyr)
df0 <- bigreadr::fread2(
  "UKBB/ukb41181.csv",
  select = c("eid", "34-0.0", "52-0.0", "22001-0.0", "21022-0.0", "189-0.0",
             paste0("22009-0.", 1:16)),
  col.names = c("eid", "year", "month", "sex", "age", "deprivation_index",
                paste0("PC", 1:16))
) %>%
  mutate(date = (year - 1900) + (month - 0.5) / 12, year = NULL, month = NULL)

ind.val <- ukb$fam$id_csv
covar <- as.matrix(df0[ind.val, -1])

ALL_METHODS <- c("lassosum2", "LDpred2", "LDpred2-low-h2",
                 "LDpred2-auto", "LDpred2-auto-rob")
ALL_CORR <- c(
  corr_UK             = "data/corr/chr",
  corr_UK_with_blocks = "data/corr/adj_with_blocks_chr",
  corr_FIN             = "data/corr_FIN/chr",
  corr_FIN_with_blocks = "data/corr_FIN/adj_with_blocks_chr",
  corr_1000G             = "data/corr_1000G_EUR/chr",
  corr_1000G_with_blocks = "data/corr_1000G_EUR/adj_with_blocks_chr"
)

bigassertr::assert_dir("results-final_FIN")

library(dplyr)
grid <- tibble::tribble(
  ~pheno,   ~phecode,
   "t1d",    "250.1",
   "t2d",    "250.2",
  "brca",    "174.1",
   "cad",    "411.4",
  "prca",      "185",
) %>%
  tidyr::expand_grid(name_corr = names(ALL_CORR), method = ALL_METHODS) %>%
  print()


library(future.apply)
plan("multisession", workers = NCORES)

grid$effects <- furrr::future_pmap(grid[1:4], function(pheno, phecode, name_corr, method) {

  # pheno <- "t2d"
  # phecode <- "250.2"
  # name_corr <- "corr_1000G"
  # method <- "lassosum2"

  res_file <- paste0("results_FIN/", pheno, "_", name_corr, "_", method, ".rds")
  if (!file.exists(res_file)) {
    print(glue::glue("/!\\ '{res_file}' IS MISSING /!\\"))
    return(rep(NA, ncol(G)))
  }

  res_file2 <- paste0("results-final_FIN/", basename(res_file))

  runonce::save_run({

    res <- readRDS(res_file)
    prs_effects <- rep(0, ncol(G))

    num_id <- readRDS(paste0("data/sumstats_FIN/", pheno, ".rds"))[["_NUM_ID_"]]

    prs_effects[num_id] <- if (grepl("LDpred2-auto", method, fixed = TRUE)) {
      all_h2 <- sapply(res, function(auto) auto$h2_est)
      h2 <- median(all_h2)
      keep <- between(all_h2, 0.7 * h2, 1.4 * h2)
      all_p <- sapply(res, function(auto) auto$p_est)
      p <- median(all_p[keep])
      keep <- keep & between(all_p, 0.5 * p, 2 * p)
      if (sum(keep) > 0) {
        rowMeans(as.matrix(sapply(res[keep], function(auto) auto$beta_est)))
      } else rep(0, length(num_id))
    } else {
      y <- readRDS("data/all_phecodes.rds")[[phecode]]
      nona <- which(!is.na(y[ind.val]) & complete.cases(covar))
      covar.nona <- covar[nona, ]
      y.val.nona <- y[ind.val[nona]]

      pred <- big_prodMat(G, res, ind.row = nona, ind.col = num_id)
      val_scores <- apply(pred, 2, function(x) {
        if (all(is.na(x))) return(NA)
        mod <- `if`(is.logical(y),
                    glm(y.val.nona ~ x + covar.nona, family = "binomial"),
                    lm(y.val.nona ~ x + covar.nona))
        summary(mod)$coef["x", 3]
      })
      res[, which.max(val_scores)]
    }

    prs_effects

  }, file = res_file2, timing = FALSE)
})

plan("sequential")


#### Evaluate final models ####

library(bigsnpr)
ukb_test <- snp_attach("data/UKBB_HM3_test.rds")
G_test <- ukb_test$genotypes
ind.test <- ukb_test$fam$id_csv
covar.test <- as.matrix(df0[ind.test, -1])

bigparallelr::set_blas_ncores(NCORES)
all_pred <- big_prodMat(G_test, do.call("cbind", grid$effects))
grid$pred <- as.list(as.data.frame(all_pred))

all_phecode <- readRDS("data/all_phecodes.rds")[unique(grid$phecode)]
plan("multisession", workers = NCORES)
grid$pcor <- furrr::future_pmap(grid[c("phecode", "pred")], function(phecode, pred) {
  if (all(is.na(pred))) return(rep(NA, 3))
  y <- all_phecode[[phecode]]
  pcor(pred, y[ind.test], covar.test)
})

plan("sequential")

ALL_POP <- c(
  corr_UK             = "UK (10,000)",
  corr_UK_with_blocks = "UK (10,000)",
  corr_FIN             = "Finnish-like (503)",
  corr_FIN_with_blocks = "Finnish-like (503)",
  corr_1000G             = "1000G European (503)",
  corr_1000G_with_blocks = "1000G European (503)"
)

ALL_PHENO <- c(
  brca = "Breast cancer (BrCa)",
  prca = "Prostate cancer (PrCa)",
  t2d = "Type 2 diabetes (T2D)",
  t1d = "Type 1 diabetes (T1D)",
  cad = "Coronary arteary disease (CAD)"
)

all_res0 <- grid %>%
  mutate(use_blocks = grepl("with_blocks", name_corr),
         pop = ALL_POP[name_corr],
         pheno = ALL_PHENO[pheno],
         method = factor(method, levels = ALL_METHODS)) %>%
  select(-effects, -pred, -phecode, -name_corr) %>%
  tidyr::unnest_wider("pcor", names_sep = "_") %>%
  mutate(across(starts_with("pcor_"), function(x) sign(x) * x^2)) %>%
  print()
# saveRDS(all_res0, "results-final_FIN/all_res.rds")



library(ggplot2)
ggplot(filter(all_res0, !use_blocks),
       aes(method, pcor_1, fill = pop)) +
  facet_wrap(~ pheno, nrow = 3, scales = "free_y") +
  bigstatsr::theme_bigstatsr() +
  scale_fill_manual(values =  c("#E69F00", "#56B4E9", "#009E73")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = pcor_2, ymax = pcor_3),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  labs(x = "Method", y = "Partial phenotypic variance explained by PGS",
       fill = "LD reference population") +
  theme(legend.position = c(0.75, 0.15),
        axis.text = element_text(size = rel(1.2 * 0.8))) +
  geom_col(data = filter(all_res0, use_blocks),
           position = position_dodge(), color = "red",
           alpha = 0, show.legend = FALSE)

# ggsave("figures/res-FIN.pdf", width = 14, height = 11)
