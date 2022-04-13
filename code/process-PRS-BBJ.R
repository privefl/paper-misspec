library(bigsnpr)
ukb <- snp_attach("data/EAS_HM3.rds")
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

id_csv <- ukb$fam$id_csv
ind.val <- which(ukb$fam$set == "val")
covar.val <- as.matrix(df0[id_csv[ind.val], -1])

ALL_METHODS <- c("lassosum2", "LDpred2", "LDpred2-low-h2",
                 "LDpred2-auto", "LDpred2-auto-rob")

ALL_CORR <- c(
  corr_JPN               = "data/corr_JPN/chr",
  corr_JPN_with_blocks   = "data/corr_JPN/adj_with_blocks_chr",
  corr_1000G             = "data/corr_1000G_EAS/chr",
  corr_1000G_with_blocks = "data/corr_1000G_EAS/adj_with_blocks_chr",
  corr_UKBB              = "data/corr_UKBB_EAS/chr",
  corr_UKBB_with_blocks  = "data/corr_UKBB_EAS/adj_with_blocks_chr"
)

bigassertr::assert_dir("results-final_BBJ")

library(dplyr)
grid <- tibble(pheno = c("height", "systolic_bp", "hdl_cholesterol", "bmi")) %>%
  tidyr::expand_grid(name_corr = names(ALL_CORR), method = ALL_METHODS) %>%
  print()


library(future.apply)
plan("multisession", workers = NCORES)

grid$effects <- furrr::future_pmap(grid[1:3], function(pheno, name_corr, method) {

  # pheno <- "height"
  # name_corr <- "corr_1000G"
  # method <- "lassosum2"

  res_file <- paste0("results_BBJ/", pheno, "_", name_corr, "_", method, ".rds")
  if (!file.exists(res_file)) {
    print(glue::glue("/!\\ '{res_file}' IS MISSING /!\\"))
    return(rep(NA, ncol(G)))
  }

  res_file2 <- paste0("results-final_BBJ/", basename(res_file))

  runonce::save_run({

    res <- readRDS(res_file)
    prs_effects <- rep(0, ncol(G))

    num_id <- readRDS(paste0("data/sumstats_BBJ/", pheno, ".rds"))[["_NUM_ID_"]]

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
      y <- readRDS("data/all_phecodes.rds")[id_csv, pheno]
      nona <- which(!is.na(y[ind.val]) & complete.cases(covar.val))
      covar.val.nona <- covar.val[nona, ]
      y.val.nona <- y[ind.val[nona]]

      pred <- big_prodMat(G, res, ind.row = ind.val[nona], ind.col = num_id)
      val_scores <- apply(pred, 2, function(x) {
        if (all(is.na(x))) return(NA)
        mod <- `if`(is.logical(y),
                    glm(y.val.nona ~ x + covar.val.nona, family = "binomial"),
                    lm(y.val.nona ~ x + covar.val.nona))
        summary(mod)$coef["x", 3]
      })
      res[, which.max(val_scores)]
    }

    prs_effects

  }, file = res_file2, timing = FALSE)
})

plan("sequential")


#### Evaluate final models ####

ind.test <- which(ukb$fam$set == "test")
covar.test <- as.matrix(df0[id_csv[ind.test], -1])

bigparallelr::set_blas_ncores(NCORES)
all_pred <- big_prodMat(G, do.call("cbind", grid$effects), ind.row = ind.test)
grid$pred <- as.list(as.data.frame(all_pred))

all_phecode <- readRDS("data/all_phecodes.rds")[id_csv, ]
grid$pcor <- purrr::pmap(grid[c("pheno", "pred")], function(pheno, pred) {
  if (all(is.na(pred))) return(rep(NA, 3))
  y <- all_phecode[[pheno]]
  pcor(pred, y[ind.test], covar.test)
})

ALL_POP <- c(
  corr_JPN               = "Japanese-like (504)",
  corr_JPN_with_blocks   = "Japanese-like (504)",
  corr_1000G             = "1000G East Asian (504)",
  corr_1000G_with_blocks = "1000G East Asian (504)",
  corr_UKBB              = "UKBB East Asian (2041)",
  corr_UKBB_with_blocks  = "UKBB East Asian (2041)"
)

ALL_PHENO <- c(
  height = "Height",
  bmi = "Body Mass Index",
  systolic_bp = "Systolic Blood Pressure",
  hdl_cholesterol = "HDL Cholesterol"
)

all_res0 <- grid %>%
  mutate(use_blocks = grepl("with_blocks", name_corr),
         pop = ALL_POP[name_corr],
         pheno = ALL_PHENO[pheno],
         method = factor(method, levels = ALL_METHODS)) %>%
  select(-effects, -pred, -name_corr) %>%
  tidyr::unnest_wider("pcor", names_sep = "_") %>%
  mutate(across(starts_with("pcor_"), function(x) sign(x) * x^2)) %>%
  print()
# saveRDS(all_res0, "results-final_BBJ/all_res.rds")


library(ggplot2)
ggplot(filter(all_res0, !use_blocks),
       aes(method, pcor_1, fill = pop)) +
  facet_wrap(~ pheno, scales = "free_y") +
  bigstatsr::theme_bigstatsr() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.02)) +
  scale_fill_manual(values =  c("#E69F00", "#56B4E9", "#009E73")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = pcor_2, ymax = pcor_3),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  labs(x = "Method", y = "Partial phenotypic variance explained by PGS",
       fill = "LD reference population") +
  theme(legend.position = "top",
        legend.text = element_text(margin = margin(r = 10, unit = "pt")),
        legend.title = element_text(margin = margin(r = 10, unit = "pt"))) +
  geom_col(data = filter(all_res0, use_blocks),
           position = position_dodge(), color = "red",
           alpha = 0, show.legend = FALSE)

# ggsave("figures/res-BBJ.pdf", width = 14, height = 10)
