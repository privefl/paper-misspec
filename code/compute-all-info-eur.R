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

set.seed(1); sampled <- slice_sample(merged, n = 50e3)

library(bigsnpr)
map_subset <- snp_attach(snp_readBGEN(
  bgenfiles   = "UKBB/bgen/ukb_imp_chr22_v3.bgen",
  list_snp_id = list(with(sampled, paste(chromosome, physical.pos, allele1, allele2, sep = "_"))),
  backingfile = tempfile(),
  ncores      = nb_cores()
))$map

library(ggplot2)
plot_grid(
  ggplot(aes(pmin(map_subset$freq, (1 - map_subset$freq)), sampled$V6), data = NULL) +
    theme_bw(14) +
    geom_point(alpha = 0.1) +
    labs(x = "MAF recomputed using all UKBB individuals", y = "MAF reported in MFI files") +
    coord_equal() +
    geom_abline(color = "red") +
    scale_x_log10(breaks = c(0.1, 0.001, 1e-5, 1e-7)) +
    scale_y_log10(breaks = c(0.1, 0.001, 1e-5, 1e-7)),

  ggplot(aes(map_subset$info, sampled$V8), data = NULL) +
    theme_bw(14) +
    geom_point(alpha = 0.2) +
    labs(x = "INFO recomputed using all UKBB individuals", y = "INFO reported in MFI files",
         color = "MAF_reported / MAF_recomputed") +
    coord_equal() +
    geom_abline(color = "red") +
    scale_color_viridis_c(trans = "log", breaks = 10^(-3:4)) +
    theme(legend.position = "top") +
    guides(color = guide_colorbar(barwidth = 12, ticks.linewidth = 2)),

  ggplot(sampled, aes(pmin(freq, (1 - freq)), V6)) +
    theme_bw(14) +
    geom_point(alpha = 0.1) +
    labs(x = "MAF recomputed in European subset", y = "MAF reported in MFI files") +
    coord_equal() +
    geom_abline(color = "red") +
    scale_x_log10(breaks = c(0.1, 0.001, 1e-5, 1e-7)) +
    scale_y_log10(breaks = c(0.1, 0.001, 1e-5, 1e-7)),

  ggplot(sampled, aes(info, V8, color = V6 / pmin(freq, (1 - freq)))) +
    theme_bw(14) +
    geom_point(alpha = 0.2) +
    labs(x = "INFO recomputed in European subset", y = "INFO reported in MFI files",
         color = expression(frac(MAF_reported, MAF_recomputed))) +
    coord_equal() +
    geom_abline(color = "red") +
    scale_color_viridis_c(trans = "log", breaks = 10^(-3:4)) +
    theme(legend.position = c(0.75, 0.15), legend.title.align = 0.5) +
    guides(color = guide_colorbar(barwidth = 12, ticks.linewidth = 1.5,
                                  direction = "horizontal", title.position = "top")),

  scale = 0.98
)
# ggsave("figures/compare-all-maf-info-chr22.png", width = 13, height = 12)


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
