library(dplyr)
library(bigreadr)

#### Get East Asian individuals close to "Japanese" individuals from the UKBB ####

df_UKBB <- runonce::save_run({

  # Country of birth
  code_country   <- filter(fread2("UKBB/coding89.tsv"), selectable == "Y")

  eid_imputed <- fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")$ID_1[-1]

  fread2(
    "UKBB/ukb41181.csv",
    select = c("eid", "20115-0.0", paste0("22009-0.", 1:16)),
    col.names = c("eid", "country", paste0("PC", 1:16))
  ) %>%
    filter(eid %in% eid_imputed) %>%
    mutate(country = factor(country, levels = code_country$coding,
                            labels = code_country$meaning))

}, file = "tmp-data/info-UKBB.rds")

center_JPN <- df_UKBB %>%
  filter(country == "Japan") %>%
  select(PC1:PC16) %>%
  bigutilsr::geometric_median()

ldist_to_JPN <- df_UKBB %>%
  select(PC1:PC16) %>%
  { sweep(., 2, center_JPN, '-')^2 } %>%
  rowSums() %>%
  sqrt() %>%
  log()

hist(ldist_to_JPN, "FD")

print_table <- function(x, min_count = 5) {
  tab <- sort(table(x, exclude = NULL), decreasing = TRUE)
  tab <- tab[tab >= min_count]
  cat(paste0(names(tab), ": ", tab), sep = " || ")
}

# Remove withdrawn and related individuals
eid_withdrawn <- scan("UKBB/w58024_20220222.csv")
rel <- fread2("UKBB/ukb58024_rel_s488264.dat")
ldist_to_JPN[df_UKBB$eid %in% c(eid_withdrawn, filter(rel, Kinship > 2^-2.5)$ID2)] <- Inf

hist(ldist_to_JPN)
length(ind <- which(ldist_to_JPN < 4.5))  # 2041
print_table(df_UKBB$country[ind])
# Hong Kong: 450 || China: 373 || Malaysia: 284 || Japan: 241 || NA: 130 || Nepal: 116 ||
# Singapore: 91 || Vietnam: 67 || Thailand: 55 || Mauritius: 34 || Myanmar (Burma): 27 ||
# South Korea: 26 || Taiwan: 25 || Indonesia: 21 || India: 12 || Brunei: 9 ||
# Caribbean: 8 || Philippines: 8 || USA: 7 || Macau (Macao): 6 || Mongolia: 6 ||
# North Korea: 6 || Brazil: 6 || South Africa: 5 || The Guianas: 5

#### Read in ####

library(bigsnpr)
sample <- fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
snp_readBGEN(
  bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
  list_snp_id = readRDS("data/list_snp_id_HM3.rds"),
  backingfile = "tmp-data/EAS_UKBB_HM3",
  ind_row     = match(df_UKBB$eid[ind], sample$ID_2),
  ncores      = nb_cores()
)

obj.bigsnp <- snp_attach("tmp-data/EAS_UKBB_HM3.rds")
G <- obj.bigsnp$genotypes
eid_csv <- bigreadr::fread2("UKBB/ukb41181.csv", select = "eid")[[1]]
set.seed(1); obj.bigsnp$fam <- data.frame(
  set    = sample(c("val", "test"), nrow(G), replace = TRUE),
  id_csv = match(df_UKBB$eid[ind], eid_csv)
)
snp_save(obj.bigsnp)

G2 <- snp_attach("tmp-data/JPN_HM3.rds")$genotypes
af_EAS <- big_scale()(G, ncores = nb_cores())$center
af_JPN  <- big_scale()(G2, ncores = nb_cores())$center
plot(af_EAS, af_JPN)  # ok


#### Prepare correlations ####

CHR <- as.integer(obj.bigsnp$map$chromosome)

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = 15, mem = "60g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("data/corr_UKBB_EAS")
bigassertr::assert_dir("data/corr_UKBB_EAS/ukbb_for_prscs")

all_final_grp <- furrr::future_map_dfr(1:22, function(chr) {

  ind.chr <- which(CHR == chr)
  MAF <- snp_MAF(G, ind.col = ind.chr, ncores = nb_cores())
  keep <- (MAF > 0.01)
  ind.chr2 <- ind.chr[keep]

  library(bigsnpr)
  corr0 <- runonce::save_run({
    POS2 <- snp_attach("data/UKBB_HM3_val.rds")$map$genetic.dist
    snp_cor(G, ind.col = ind.chr2, infos.pos = POS2[ind.chr2],
            size = 3 / 1000, ncores = nb_cores())
  }, file = paste0("data/corr_UKBB_EAS/chr", chr, ".rds"))


  # find nearly independent LD blocks
  m <- length(ind.chr2)
  (SEQ <- round(seq_log(m / 30, m / 5, length.out = 20)))
  splits <- snp_ldsplit(corr0, thr_r2 = 0.03, min_size = 50, max_size = SEQ, max_r2 = 0.15)
  splits$cost2 <- sapply(splits$all_size, function(sizes) sum(sizes^2))

  library(dplyr)
  best_split <- splits %>%
    arrange(cost2 * sqrt(5 + cost)) %>%
    print() %>%
    slice(1) %>%
    print()

  (all_size <- best_split$all_size[[1]])
  best_grp <- rep(seq_along(all_size), all_size)

  runonce::save_run({
    corr0T <- as(corr0, "dgTMatrix")
    corr0T@x <- ifelse(best_grp[corr0T@i + 1L] == best_grp[corr0T@j + 1L], corr0T@x, 0)
    as(Matrix::drop0(corr0T), "symmetricMatrix")
  }, file = paste0("data/corr_UKBB_EAS/adj_with_blocks_chr", chr, ".rds"))


  # LD for PRS-CS (based on the same data and LD blocks)
  hdf5_file <- paste0("data/corr_UKBB_EAS/ukbb_for_prscs/ldblk_ukbb_chr", chr, ".hdf5")
  runonce::skip_run_if({
    ind_block <- split(ind.chr2, best_grp)
    bigparallelr::set_blas_ncores(nb_cores())
    for (ic in seq_along(ind_block)) {
      ind <- ind_block[[ic]]
      ld <- big_cor(G, ind.col = ind)[]
      rhdf5::h5write(list(ldblk = ld, snplist = paste0("snp", ind)),
                     file = hdf5_file, name = paste0("blk_", ic), level = 9)
    }
  }, files = hdf5_file)


  # return
  tibble(best_split, ind = list(ind.chr2))
})


# verif
plot(all_final_grp$n_block)
plot(all_final_grp$cost)
plot(all_final_grp$cost2)
# saveRDS(all_final_grp, "data/corr_UKBB_EAS/all_final_grp.rds")

sum(file.size(paste0("data/corr_UKBB_EAS/adj_with_blocks_chr", 1:22, ".rds"))) /
  sum(file.size(paste0("data/corr_UKBB_EAS/chr", 1:22, ".rds")))
# 48.5%


# need also some .bim and some snpinfo for PRS-CS
ind_keep <- unlist(all_final_grp$ind)
snpinfo <- transmute(obj.bigsnp$map, chromosome, marker.ID = paste0("snp", row_number()),
                     genetic.dist = 0, physical.pos, allele1 = "A", allele2 = "C")[ind_keep, ]
bigsnpr:::write.table2(snpinfo, "data/corr_UKBB_EAS/ukbb_for_prscs/for_prscs.bim")

snpinfo2 <- transmute(snpinfo, CHR = chromosome, SNP = marker.ID,
                      BP = physical.pos, A1 = allele1, A2 = allele2)
snpinfo2$MAF <- snp_MAF(G, ind.col = ind_keep, ncores = nb_cores())
bigreadr::fwrite2(snpinfo2, "data/corr_UKBB_EAS/ukbb_for_prscs/snpinfo_ukbb_hm3", sep = "\t")
