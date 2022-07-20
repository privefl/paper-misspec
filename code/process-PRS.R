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
                 "LDpred2-auto", "LDpred2-auto-rob", "PRS-CS-auto")

bigassertr::assert_dir("results-final")
bigassertr::assert_dir("results-final/sumstats")

library(dplyr)
grid <- tibble::tribble(
  ~pheno,                               ~phecode,
  paste0("t1d_", c("affy", "illu")),    "250.1",
  paste0("brca_", c("onco", "icogs")),  "174.1",
  "cad",                                "411.4",
  "prca",                               "185",
  "mdd",                                "296.22",
) %>%
  tidyr::unnest("pheno") %>%
  mutate(whichN = NA) %>%
  tidyr::expand_grid(use_info = c(FALSE, TRUE)) %>%
  add_row(pheno = "vitaminD", phecode = "vitaminD", whichN = c("trueN", "maxN"), use_info = FALSE) %>%
  tidyr::expand_grid(
    qc = c("noqc", "qc1", "qc2"),
    name_corr = c("data/corr/chr", "data/corr/adj_with_blocks_chr"),
    method = ALL_METHODS) %>%
  mutate(name_corr = ifelse(method == "PRS-CS-auto", NA_character_, name_corr)) %>%
  unique() %>%
  print()
# 528 different PGS models


library(future.apply)
plan("multisession", workers = NCORES)

grid$effects <- furrr::future_pmap(grid[1:7], function(
  pheno, phecode, qc, use_info, whichN, name_corr, method) {

  # pheno <- "t1d_illu"
  # phecode <- "250.1"
  # qc <- "qc1"
  # use_info <- TRUE
  # name_corr <- "data/corr/adj_with_blocks_chr"
  # method <- "PRS-CS-auto"

  basename <- paste0(pheno, "_", qc)
  sumstats_file <- paste0("data/sumstats/", basename, ".rds")

  if (method != "PRS-CS-auto") {

    suffix_N <- `if`(is.na(whichN), "", paste0("_", whichN))
    basename <- paste0("results/sumstats/", basename, suffix_N)
    if (grepl("blocks", name_corr, fixed = TRUE))
      basename <- paste0(basename, "_adj_with_blocks")

    res_file <- paste0(basename, "_", method, ".rds")
    if (!file.exists(res_file)) {
      print(glue::glue("/!\\ '{res_file}' IS MISSING /!\\"))
      return(rep(NA, ncol(G)))
    }
    res <- readRDS(res_file)
  } else {
    res_files <- paste0("results_prscs/sumstats/", pheno, "_", qc, "_chr", 1:22, ".rds")
    for (res_file in res_files) {
      if (!file.exists(res_file)) {
        print(glue::glue("/!\\ '{res_file}' IS MISSING /!\\"))
        return(rep(NA, ncol(G)))
      }
    }
    res <- unlist(lapply(res_files, function(res_file) readRDS(res_file)$V6))
  }

  res_file2 <- paste0("results-final/sumstats/", basename(basename), "_", method,
                      `if`(use_info, "_info", ""), ".rds")

  runonce::save_run({

    num_id <- readRDS(sumstats_file)[["_NUM_ID_"]]
    sqrt_info <- `if`(use_info, sqrt(readRDS(sumstats_file)[["info"]]), rep(1, length(num_id)))
    if (pheno == "vitaminD") sqrt_info <- -sqrt_info

    prs_effects <- rep(0, ncol(G))
    prs_effects[num_id] <- if (grepl("LDpred2-auto", method, fixed = TRUE)) {
      range <- sapply(res, function(auto) diff(range(auto$corr_est)))
      keep <- (range > (0.9 * quantile(range, 0.9)))
      rowMeans(sapply(res[keep], function(auto) auto$beta_est)) * sqrt_info
    } else if (method == "PRS-CS-auto") {
      res * sqrt_info
    } else {
      y <- readRDS("data/all_phecodes.rds")[[phecode]]
      nona <- which(!is.na(y[ind.val]) & complete.cases(covar))
      covar.nona <- covar[nona, ]
      y.val.nona <- y[ind.val[nona]]

      res <- sweep(res, 1, sqrt_info, '*')
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

all_res0 <- grid %>%
  mutate(use_blocks = method == "PRS-CS-auto" | grepl("with_blocks", name_corr),
         method = factor(method, levels = ALL_METHODS)) %>%
  select(-effects, -pred, -phecode, -name_corr) %>%
  tidyr::unnest_wider("pcor", names_sep = "_") %>%
  mutate(across(starts_with("pcor_"), function(x) sign(x) * x^2)) %>%
  filter(method != "PRS-CS-auto" | whichN %in% c(NA, "maxN")) %>%
  print()
# saveRDS(all_res0, "results-final/sumstats/all_res.rds")


for (PHENO in setdiff(unique(all_res0$pheno), "vitaminD")) {

  all_res <- filter(all_res0, pheno == PHENO) %>%
    mutate(use_info = ifelse(use_info, "Yes", "No"))

  library(ggplot2)
  ggplot(filter(all_res, method == "PRS-CS-auto" | !use_blocks),
         aes(method, pcor_1, fill = paste(qc, use_info, sep = " - "))) +
    facet_wrap(~ pheno, ncol = 1) +
    bigstatsr::theme_bigstatsr(0.9) +
    scale_fill_manual(values =  c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")) +
    geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
    geom_errorbar(aes(ymin = pcor_2, ymax = pcor_3),
                  position = position_dodge(width = 0.9),
                  color = "black", width = 0.2, size = 1) +
    labs(x = "Method", y = "Partial phenotypic variance explained by PGS",
         fill = "Which QC? - Correction using INFO?") +
    theme(legend.position = "top",
          legend.text = element_text(margin = margin(r = 10, unit = "pt")),
          legend.title = element_text(margin = margin(r = 10, unit = "pt"))) +
    geom_col(data = filter(all_res, use_blocks),
             position = position_dodge(), color = "red",
             alpha = 0, show.legend = FALSE)

  ggsave(paste0("figures/res-", PHENO, ".pdf"), width = 9.5, height = 6.5)
}


PHENO <- "vitaminD"
all_res <- filter(all_res0, grepl(PHENO, pheno))

library(ggplot2)
ggplot(filter(all_res, method == "PRS-CS-auto" | !use_blocks),
       aes(method, pcor_1, fill = paste(qc, whichN, sep = " - "))) +
  facet_wrap(~ pheno, ncol = 1) +
  bigstatsr::theme_bigstatsr(0.9) +
  scale_fill_manual(values =  c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")) +
  geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
  geom_errorbar(aes(ymin = pcor_2, ymax = pcor_3),
                position = position_dodge(width = 0.9),
                color = "black", width = 0.2, size = 1) +
  labs(x = "Method", y = "Partial phenotypic variance explained by PGS",
       fill = "Which QC? - Which N?") +
  theme(legend.position = "top",
        legend.text = element_text(margin = margin(r = 10, unit = "pt")),
        legend.title = element_text(margin = margin(r = 10, unit = "pt"))) +
  geom_col(data = filter(all_res, use_blocks),
           position = position_dodge(), color = "red",
           alpha = 0, show.legend = FALSE)

ggsave("figures/res-vitaminD.pdf", width = 9.5, height = 6.5)
