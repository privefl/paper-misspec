library(bigsnpr)
ukb <- snp_attach("data/UKBB_HM3_val.rds")
G <- ukb$genotypes
CHR <- ukb$map$chromosome
POS2 <- ukb$map$genetic.dist
INFO <- ukb$map$info

(NCORES <- nb_cores())

chr <- 22
ind.chr <- which(CHR == chr)

corr <- runonce::save_run(
  snp_cor(G, ind.col = ind.chr, infos.pos = POS2[ind.chr],
          size = 3 / 1000, ncores = NCORES),
  file = paste0("data/corr/chr", chr, ".rds")
)

corr2 <- runonce::save_run(
  snp_cor(G, ind.col = ind.chr, infos.pos = POS2[ind.chr],
          size = 3 / 1000, info = INFO[ind.chr], ncores = NCORES),
  file = paste0("data/corr/adj_chr", chr, ".rds")
)

ind_csv <- ukb$fam$id_csv
all_id <- bigreadr::fread2("UKBB/ukb41181.csv", select = "eid")[[1]]
sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
ind.indiv <- match(all_id[ind_csv], sample$ID_2)

list_snp_id <- with(ukb$map[ind.chr, ],
                    split(paste(chromosome, physical.pos, allele1, allele2, sep = "_"),
                          as.integer(chromosome)))

for (ic in 1:10) {

  print(ic)

  rds <- snp_readBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen"),
    list_snp_id = list_snp_id,
    backingfile = tempfile(tmpdir = "tmp-data/sfbm/"),
    ind_row     = ind.indiv,
    read_as     = "random"
  )

  if (ic == 1) {
    corr3 <- snp_cor(snp_attach(rds)$genotypes, infos.pos = POS2[ind.chr],
                     size = 3 / 1000, ncores = NCORES)
  } else {
    corr3 <- corr3 + snp_cor(snp_attach(rds)$genotypes, infos.pos = POS2[ind.chr],
                             size = 3 / 1000, ncores = NCORES)
  }

  file.remove(c(rds, sub("\\.rds$", ".bk", rds)))
}

corr3 <- corr3 / ic


sub <- sample(length(ind.chr), 1000, prob = 1 / INFO[ind.chr]^100)
hist(INFO[ind.chr][sub])

ind <- Matrix::which(corr[sub, sub] != 0, arr.ind = TRUE)
ind2 <- cbind(sub[ind[, 1]], sub[ind[, 2]])
INFO_ind2 <- sqrt(INFO[ind.chr[ind2[, 1]]] * INFO[ind.chr[ind2[, 2]]])

library(ggplot2)
source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
plot_grid2(list(
  qplot(corr3[ind2], corr[ind2], color = INFO_ind2, alpha = I(0.5)) +
    scale_color_viridis_c(direction = -1) +
    theme_bw() +
    labs(x = "Correlations from multiple imputation",
         y = "Correlations from imputed data",
         color = "INFO") +
    guides(color = guide_colorbar(barheight = 15, ticks.linewidth = 2)) +
    geom_abline(col = "red"),
  qplot(corr3[ind2], corr2[ind2], color = INFO_ind2, alpha = I(0.5)) +
    scale_color_viridis_c(direction = -1) +
    # coord_equal() +
    theme_bw() +
    labs(x = "Correlations from multiple imputation",
         y = "Correlations from imputed data (with correction)") +
    geom_abline(col = "red")
), scale = 0.95, title_ratio = 0, legend_ratio = 0.1)
# ggsave("figures/compare-cor.png", width = 10, height = 5)
