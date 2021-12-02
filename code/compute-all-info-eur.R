## Get all SNP IDs and split them in blocks
library(dplyr)

all_snp_id <- runonce::save_run(
  purrr::map_dfr(1:22, function(chr) {

    print(chr)

    bgifile <- paste0("UKBB/bgen/ukb_imp_chr", chr, "_v3.bgen.bgi")
    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgifile)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)

    snp_id <- tbl(db_con, "Variant") %>%
      select(chr = chromosome, pos = position, a1 = allele1, a0 = allele2) %>%
      collect() %>%
      with(paste(chr, pos, a1, a0, sep = "_"))

    all_splits <- bigparallelr::split_vec(snp_id, block_len = 50e3)

    tibble(bgifile, snp_id = all_splits)
  }),
  file = "tmp-data/split_all_variants.rds"
)


## Loop over all blocks
grid <- mutate(all_snp_id, id = row_number())  # 1874 parts

library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("tmp-data/bgen")
bigassertr::assert_dir("all_info")

furrr::future_pwalk(grid, function(bgifile, snp_id, id) {

  runonce::save_run({
    tmp <- tempfile(tmpdir = "tmp-data/bgen")
    on.exit(file.remove(paste0(tmp, c(".bk", ".rds"))), add = TRUE)

    library(bigsnpr)
    map <- snp_attach(snp_readBGEN(
      bgenfiles   = sub("\\.bgi$", "", bgifile),
      list_snp_id = list(snp_id),
      backingfile = tmp,
      ind_row     = readRDS("tmp-data/eur_ind_bgen.rds"),
      ncores      = nb_cores()
    ))$map
  }, file = paste0("all_info/part", id, ".rds"))
})


## Look at chr22
all_info <- purrr::map_dfr(grep("chr22", grid$bgifile), function(id) {
  readRDS(paste0("all_info/part", id, ".rds"))
})

mfi <- bigreadr::fread2("UKBB/mfi/ukb_mfi_chr22_v3.txt")
merged <- inner_join(all_info, mfi,
                     by = c(rsid = "V2", allele1 = "V4", allele2 = "V5"))

sampled <- slice_sample(merged, n = 50e3)
library(ggplot2)
ggplot(sampled, aes(pmin(freq, (1 - freq)), V6)) +
  theme_bw(14) +
  geom_point(alpha = 0.1) +
  labs(x = "MAF recomputed", y = "MAF reported") +
  coord_equal() +
  geom_abline(color = "red") +
  scale_x_log10() + scale_y_log10()
# ggsave("figures/compare-all-maf-chr22.png", width = 8, height = 7)

ggplot(sampled, aes(info, V8, color = V6 / pmin(freq, (1 - freq)))) +
  theme_bw(14) +
  geom_point(alpha = 0.2) +
  labs(x = "INFO recomputed", y = "INFO reported",
       color = "MAF_reported / MAF_recomputed") +
  coord_equal() +
  geom_abline(color = "red") +
  scale_color_viridis_c(trans = "log", breaks = 10^(-3:4)) +
  theme(legend.position = "top") +
  guides(color = guide_colorbar(barwidth = 12, ticks.linewidth = 2))
# ggsave("figures/compare-all-info-chr22.png", width = 7.5, height = 8)


## Write CSVs to export
bigassertr::assert_dir("export_info")

for (chr in 1:22) {

  print(chr)

  bgifile <- paste0("UKBB/bgen/ukb_imp_chr", chr, "_v3.bgen.bgi")
  parts <- which(all_snp_id$bgifile == bgifile)

  all_info <- purrr::map_dfr(paste0("all_info/part", parts, ".rds"), readRDS)
  all_info$chromosome <- as.integer(all_info$chromosome)

  bigreadr::fwrite2(
    all_info,
    file = paste0("export_info/ukbb_eur_info_chr", chr, ".csv.gz"),
    compress = "gzip")
}
