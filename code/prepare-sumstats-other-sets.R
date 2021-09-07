# dput(list.files("code/prepare-sumstats/"))
files <- c(
  "BrCa-icogs.R",
  "BrCa-onco.R",
  "CAD.R",
  "MDD.R",
  "PrCa.R",
  "T1D-affy.R",
  "T1D-illu.R"
)


bigassertr::assert_dir("data/sumstats-maxtag")

for (file in files) {

  rm(list = setdiff(ls(), c("file", "files")))

  code <- readLines(file.path("code/prepare-sumstats", file))

  library(magrittr)
  new_code <- code %>%
    sub('"data/UKBB_HM3_val.rds"', '"data/UKBB_maxtag_val.rds"', ., fixed = TRUE) %>%
    sub('"data/sumstats/', '"data/sumstats-maxtag/', ., fixed = TRUE) %>%
    sub('# saveRDS(filter(info_snp2', 'saveRDS(filter(info_snp2', ., fixed = TRUE) %>%
    print()

  stopifnot(sum(new_code != code) == 3)

  eval(parse(text = new_code))
}


bigassertr::assert_dir("data/sumstats-clump")

map <- readRDS("data/UKBB_large_val.rds")$map
ind_clump <- readRDS("data/ind_keep_large_clump.rds")
saveRDS(list(map = map[ind_clump, ]), "tmp-data/fake_for_clump.rds")

for (file in files) {

  rm(list = setdiff(ls(), c("file", "files")))

  code <- readLines(file.path("code/prepare-sumstats", file))

  library(magrittr)
  new_code <- code %>%
    sub('snp_attach(', 'readRDS(', ., fixed = TRUE) %>%
    sub('"data/UKBB_HM3_val.rds"', '"tmp-data/fake_for_clump.rds"', ., fixed = TRUE) %>%
    sub('"data/sumstats/', '"data/sumstats-clump/', ., fixed = TRUE) %>%
    sub('# saveRDS(filter(info_snp2', 'saveRDS(filter(info_snp2', ., fixed = TRUE) %>%
    print()

  stopifnot(sum(new_code != code) == 3)

  eval(parse(text = new_code))
}
