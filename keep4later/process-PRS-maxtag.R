library(bigsnpr)
ukb <- snp_attach("data/UKBB_maxtag_val.rds")
G <- ukb$genotypes

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

bigassertr::assert_dir("results-maxtag-final")

library(dplyr)
grid <- tibble::tribble(
  ~pheno,                               ~phecode,
  paste0("t1d_", c("affy", "illu")),    "250.1",
  paste0("brca_", c("onco", "icogs")),  "174.1",
  "cad",                                "411.4",
  "prca",                               "185",
  "mdd",                                "296.22",
) %>%
  tidyr::expand_grid(
    qc = c("qc1", "qc2"),
    use_info = c(FALSE, TRUE),
    name_corr = c("data/corr-large/maxtag_adj_chr",
                  "data/corr-large/maxtag_adj_with_blocks_chr"),
    method = ALL_METHODS) %>%
  tidyr::unnest("pheno") %>%
  print()


library(future.apply)
plan("multisession", workers = 10)

grid$effects <- furrr::future_pmap(grid[1:6], function(
  pheno, phecode, qc, use_info, name_corr, method) {

  # pheno <- "t1d_illu"
  # phecode <- "250.1"
  # qc <- "qc1"
  # use_info <- TRUE
  # name_corr <- "data/corr-large/maxtag_adj_with_blocks_chr"
  # method <- "lassosum2"

  basename <- paste0(pheno, "_", qc)
  sumstats_file <- paste0("data/sumstats-maxtag/", basename, ".rds")
  if (grepl("blocks", name_corr, fixed = TRUE))
    basename <- paste0(basename, "_adj_with_blocks")
  basename <- paste0("results-maxtag/", basename)

  res_file <- paste0(basename, "_", method, ".rds")
  if (!file.exists(res_file)) {
    print(glue::glue("/!\\ '{res_file}' IS MISSING /!\\"))
    return(rep(NA, ncol(G)))
  }

  res_file2 <- paste0("results-maxtag-final/", basename(basename), "_", method,
                      `if`(use_info, "_info", ""), ".rds")

  runonce::save_run({

    res <- readRDS(res_file)
    prs_effects <- rep(0, ncol(G))

    sumstats <- readRDS(sumstats_file)[c("_NUM_ID_", "info")]

    num_id <- sumstats[[1]]
    sqrt_info <- `if`(use_info, sqrt(sumstats[[2]]), rep(1, length(num_id)))

    prs_effects[num_id] <- if (grepl("LDpred2-auto", method, fixed = TRUE)) {
      all_h2 <- sapply(res, function(auto) auto$h2_est)
      h2 <- median(all_h2)
      keep <- between(all_h2, 0.7 * h2, 1.4 * h2)
      all_p <- sapply(res, function(auto) auto$p_est)
      p <- median(all_p[keep])
      keep <- keep & between(all_p, 0.5 * p, 2 * p)
      rowMeans(sapply(res[keep], function(auto) auto$beta_est)) * sqrt_info
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
ind.test <- snp_attach("data/UKBB_HM3_test.rds")$fam$id_csv
covar.test <- as.matrix(df0[ind.test, -1])

all_id <- bigreadr::fread2("UKBB/ukb41181.csv", select = "eid")[[1]]
sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]

all_pred <- runonce::save_run(
  snp_prodBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = readRDS("data/UKBB_maxtag_snp_id.rds"),
    beta        = do.call("cbind", grid$effects),
    ind_row     = match(all_id[ind.test], sample$ID_2),
    ncores      = nb_cores()
  ),
  file = "results-maxtag-final/all_pred.rds"
) # 73 min

grid$pred <- as.list(as.data.frame(all_pred))

all_phecode <- readRDS("data/all_phecodes.rds")[unique(grid$phecode)]
plan("multisession", workers = nb_cores())
grid$pcor <- furrr::future_pmap(grid[c("phecode", "pred")], function(phecode, pred) {
  if (all(is.na(pred))) return(rep(NA, 3))
  y <- all_phecode[[phecode]]
  pcor(pred, y[ind.test], covar.test)
})

plan("sequential")

all_res0 <- grid %>%
  mutate(use_blocks = grepl("with_blocks", name_corr),
         method = factor(method, levels = ALL_METHODS)) %>%
  select(-effects, -pred, -phecode, -name_corr) %>%
  tidyr::unnest_wider("pcor", names_sep = "_") %>%
  print()
# saveRDS(all_res0, "results-maxtag-final/all_res.rds")


# for (PHENO in c("t1d", "brca", "mdd", "prca", "cad")) {
#
#   all_res <- filter(all_res0, grepl(PHENO, pheno))
#
#   library(ggplot2)
#   ggplot(filter(all_res, !use_blocks),
#          aes(method, pcor_1, fill = paste(qc, use_info))) +
#     facet_wrap(~ pheno, nrow = 1) +
#     bigstatsr::theme_bigstatsr(0.8) +
#     scale_fill_manual(values =  c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
#     geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
#     geom_errorbar(aes(ymin = pcor_2, ymax = pcor_3),
#                   position = position_dodge(width = 0.9),
#                   color = "black", width = 0.2, size = 1) +
#     labs(x = "Method", y = "Partial correlation between PGS and phenotype",
#          fill = "Which QC? Correction using INFO?") +
#     theme(legend.position = "top") +
#     geom_col(data = filter(all_res, use_blocks),
#              position = position_dodge(), color = "red",
#              alpha = 0, show.legend = FALSE)
#
#   ggsave(paste0("figures/res-maxtag-", PHENO, ".pdf"),
#          width = `if`(nrow(all_res) < 50, 9, 13), height = 6.5)
# }
