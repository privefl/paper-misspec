library(bigsnpr)
ukb <- snp_attach("data/UKBB_HM3_val.rds")
G <- ukb$genotypes
CHR <- ukb$map$chromosome
POS2 <- ukb$map$genetic.dist
INFO <- ukb$map$info
(NCORES <- nb_cores())

chr <- 1
length(ind.chr <- which(CHR == chr & INFO < 0.85))  # 214


## From 190 individuals in 1000G
obj.bed <- bed(download_1000G("tmp-data"))
fam2 <- bigreadr::fread2(paste0(obj.bed$prefix, ".fam2"))
ind_row <- which(fam2$Population %in% c("CEU", "GBR"))  # only 190 individuals
matched <- snp_match(setNames(obj.bed$map[c(1, 4:6, 3)], c("chr", "pos", "a1", "a0", "beta")),
                     setNames(cbind(ukb$map[ind.chr, c(1, 4:6)], POS2[ind.chr]),
                              c("chr", "pos", "a1", "a0", "pos2")))
ind_col <- matched$`_NUM_ID_.ss`

corr_1kg <- bed_cor(obj.bed, ind.row = ind_row, ind.col = ind_col, fill.diag = FALSE,
                    infos.pos = matched$pos2, size = 3 / 1000, ncores = NCORES)


## From the UKBB dosage data
ind.chr2 <- ind.chr[matched$`_NUM_ID_`]
hist(INFO[ind.chr2])

corr_dosage <- snp_cor(G, ind.col = ind.chr2, infos.pos = POS2[ind.chr2],
                       size = 3 / 1000, ncores = NCORES)

## From multiple imputation
ind_csv <- ukb$fam$id_csv
all_id <- bigreadr::fread2("UKBB/ukb41181.csv", select = "eid")[[1]]
sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
ind.indiv <- match(all_id[ind_csv], sample$ID_2)

list_snp_id <- with(ukb$map[ind.chr2, ],
                    split(paste(chromosome, physical.pos, allele1, allele2, sep = "_"),
                          as.integer(chromosome)))

for (ic in 1:20) {

  print(ic)

  rds <- snp_readBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen"),
    list_snp_id = list_snp_id,
    backingfile = tempfile(tmpdir = "tmp-data/sfbm/"),
    ind_row     = ind.indiv,
    read_as     = "random",
    ncores      = NCORES
  )

  if (ic == 1) {
    corr_MI <- snp_cor(snp_attach(rds)$genotypes, infos.pos = POS2[ind.chr2],
                       size = 3 / 1000, ncores = NCORES)
  } else {
    corr_MI <- corr_MI + snp_cor(snp_attach(rds)$genotypes, infos.pos = POS2[ind.chr2],
                                 size = 3 / 1000, ncores = NCORES)
  }

  file.remove(c(rds, sub("\\.rds$", ".bk", rds)))
}

corr_MI <- corr_MI / ic



ind <- Matrix::which(corr_1kg != 0, arr.ind = TRUE)
INFO_ind <- sqrt(INFO[ind.chr2[ind[, 1]]] * INFO[ind.chr2[ind[, 2]]])

library(ggplot2)
source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
plot_grid2(list(
  qplot(corr_1kg[ind], corr_dosage[ind], color = INFO_ind, alpha = I(0.5)) +
    scale_color_viridis_c(direction = 1) +
    theme_bw() +
    labs(x = "Correlations from the 1000G data",
         y = "Correlations from imputed dosage data (UKBB)",
         color = "INFO") +
    guides(color = guide_colorbar(barheight = 15, ticks.linewidth = 2)) +
    geom_abline(col = "red", linetype = 2),
  qplot(corr_1kg[ind], corr_MI[ind], color = INFO_ind, alpha = I(0.5)) +
    scale_color_viridis_c(direction = 1) +
    theme_bw() +
    labs(x = "Correlations from the 1000G data",
         y = "Correlations from multiple imputation (UKBB)") +
    guides(color = guide_colorbar(barheight = 15, ticks.linewidth = 2)) +
    geom_abline(col = "red", linetype = 2)
), scale = 0.95, title_ratio = 0, legend_ratio = 0.1)
# ggsave("figures/compare-cor.png", width = 10.5, height = 5)
