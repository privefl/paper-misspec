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
POS2 <- snp_attach("data/UKBB_HM3_val.rds")$map$genetic.dist

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = 15, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("data/corr_UKBB_EAS")

all_final_grp <- furrr::future_map_dfr(1:22, function(chr) {

  ind.chr <- which(CHR == chr)

  library(bigsnpr)
  corr0 <- runonce::save_run({
    corr0 <- snp_cor(G, ind.col = ind.chr, infos.pos = POS2[ind.chr],
                     size = 3 / 1000, ncores = nb_cores())
    mean(is.na(corr0@x))  # 0.001 for chr 22
    corr0@x <- ifelse(is.na(corr0@x) | is.infinite(corr0@x), 0, corr0@x)
    corr0
  }, file = paste0("data/corr_UKBB_EAS/chr", chr, ".rds"))

  all_splits <- runonce::save_run({

    SEQ <- round(seq_log(1000 + ncol(corr0) / 50,
                         5000 + ncol(corr0) / 10,
                         length.out = 10))

    library(furrr)
    plan("multisession", workers = 5)
    options(future.globals.maxSize = 10e9)

    splits <- future_map_dfr(SEQ, function(max_size) {
      res <- snp_ldsplit(corr0, thr_r2 = 0.05, min_size = 100,
                         max_size = max_size, max_K = 200)
      res$max_size <- max_size
      res
    })
  }, file = paste0("tmp-data/split_corr_UKBB_EAS_chr", chr, ".rds"))


  library(ggplot2)
  qplot(data = all_splits, perc_kept, cost, color = as.factor(max_size)) +
    theme_bw(12) +
    theme(legend.position = "top") + guides(colour = guide_legend(nrow = 1)) +
    labs(y = "Sum of squared correlations outside blocks (log-scale)",
         x = "% of non-zero values kept", color = "Maximum block size") +
    ylim(0, min(quantile(all_splits$cost, 0.6), 2000)) +
    geom_vline(xintercept = 0.6, linetype = 3) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2))


  library(dplyr)
  final_split <- all_splits %>%
    filter(perc_kept < `if`(chr == 6, 0.7, 0.6)) %>%
    arrange(cost) %>%
    slice(1)

  final_grp <- final_split$block_num[[1]]

  corr0T <- as(corr0, "dgTMatrix")
  corr0T@x <- ifelse(final_grp[corr0T@i + 1L] == final_grp[corr0T@j + 1L], corr0T@x, 0)

  corr0T %>%
    Matrix::drop0() %>%
    as("symmetricMatrix") %>%
    saveRDS(paste0("data/corr_UKBB_EAS/adj_with_blocks_chr", chr, ".rds"))

  final_split
})

plot(all_final_grp$n_block)
plot(all_final_grp$cost)  # had to switch to perc=0.75 for chr6 based on this plots