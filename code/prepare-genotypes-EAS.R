library(dplyr)
library(bigreadr)

#### Get some East Asian individuals from the UKBB ####

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

all_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)

center_EAS <- all_centers %>%
  filter(Ancestry == "China") %>%
  select(-Ancestry) %>%
  unlist()

ldist <- df_UKBB %>%
  select(PC1:PC16) %>%
  { sweep(., 2, center_EAS, '-')^2 } %>%
  rowSums() %>%
  sqrt() %>%
  log()

hist(ldist[ldist < 5], "FD")

print_table <- function(x, min_count = 5) {
  tab <- sort(table(x, exclude = NULL), decreasing = TRUE)
  tab <- tab[tab >= min_count]
  cat(paste0(names(tab), ": ", tab), sep = " || ")
}

# Remove withdrawn and related individuals
eid_withdrawn <- scan("UKBB/w58024_20220222.csv")
rel <- fread2("UKBB/ukb58024_rel_s488264.dat")
ldist[df_UKBB$eid %in% c(eid_withdrawn, filter(rel, Kinship > 2^-2.5)$ID2)] <- Inf

ind <- which(ldist < 3.8)  # 1777 indivs
print_table(df_UKBB$country[ind])
# Hong Kong: 441 || China: 368 || Malaysia: 273 || Japan: 236 || NA: 105 || Singapore: 86 ||
# Vietnam: 43 || Mauritius: 31 || Thailand: 28 || South Korea: 26 || Taiwan: 25 ||
# Indonesia: 21 || Myanmar (Burma): 13 || Brunei: 8 || Caribbean: 7 || India: 7 || USA: 7 ||
# Macau (Macao): 6 || North Korea: 6 || Brazil: 6 || Philippines: 5 || The Guianas: 5


library(bigsnpr)
sample <- fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
snp_readBGEN(
  bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
  list_snp_id = readRDS("data/list_snp_id_HM3.rds"),
  backingfile = "data/EAS_HM3",
  ind_row     = match(df_UKBB$eid[ind], sample$ID_2),
  ncores      = nb_cores()
)

obj.bigsnp <- snp_attach("data/EAS_HM3.rds")
G <- obj.bigsnp$genotypes
eid_csv <- bigreadr::fread2("UKBB/ukb41181.csv", select = "eid")[[1]]
set.seed(1); obj.bigsnp$fam <- data.frame(
  set    = sample(c("val", "test"), nrow(G), replace = TRUE),
  id_csv = match(df_UKBB$eid[ind], eid_csv)
)
snp_save(obj.bigsnp)
