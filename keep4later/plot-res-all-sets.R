library(dplyr)

all_res0 <- bind_rows(
  bind_cols(readRDS("results-final/all_res.rds"), set = "HM3"),
  bind_cols(readRDS("results-maxtag-final/all_res.rds"), set = "maxtag"),
  bind_cols(readRDS("results-clump-final/all_res.rds"), set = "clump")
) %>%
  mutate(set = factor(set, c("HM3", "maxtag", "clump")))

params <- tribble(
  ~pheno, ~lo,  ~up,
  "prca", 0.12, 0.15,
  "brca", 0.11, 0.13,
  "cad",  0.06, 0.12,
  "mdd",  0.02, 0.04,
  "t1d",  0.04, 0.06,
)

for (PHENO in c("t1d", "brca", "mdd", "prca", "cad")) {

  lo <- filter(params, pheno == PHENO)$lo

  all_res <- filter(all_res0, grepl(PHENO, pheno))

  library(ggplot2)
  ggplot(filter(all_res, !use_blocks),
         aes(method, pcor_1, fill = set)) +
    scale_y_continuous(breaks = seq(0, 0.2, by = 0.02),
                       minor_breaks = seq(0, 0.2, by = 0.01),
                       limits = c(lo, NA), oob = scales::rescale_none) +
    facet_grid(qc + use_info ~ pheno) +
    bigstatsr::theme_bigstatsr(0.8) +
    scale_fill_manual(values =  c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
    geom_col(position = position_dodge(), alpha = 0.6, color = "black", size = 1) +
    geom_errorbar(aes(ymin = pcor_2, ymax = pcor_3),
                  position = position_dodge(width = 0.9),
                  color = "black", width = 0.2, size = 1) +
    labs(x = "Method", y = "Partial correlation between PGS and phenotype",
         fill = "Which set of variants?") +
    theme(legend.position = "top") +
    geom_col(data = filter(all_res, use_blocks),
             position = position_dodge(), color = "red",
             alpha = 0, show.legend = FALSE)

  ggsave(paste0("figures/res-all-", PHENO, ".pdf"),
         width = `if`(nrow(all_res) < 50, 9, 13), height = 6.5)
}
