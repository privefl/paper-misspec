#### Code from https://doi.org/10.1101/2021.02.05.21251061 ####

library(bigreadr)
library(dplyr)
csv <- "UKBB/ukb41181.csv"
sex <- fread2(csv, select = "22001-0.0")[[1]]

df_ICD10 <- fread2(csv, colClasses = "character", select = c(
  paste0("40001-", 0:1, ".0"),     # death
  paste0("40002-0.", 0:13),        # death
  paste0("40002-1.", 0:13),        # death
  paste0("40006-", 0:16, ".0"),    # cancer
  paste0("41270-0.", 0:212),       # diagnosis
  paste0("41201-0.", 0:21)         # external cause
))

# Non-cancer illness code, self-reported
# coding609 <- fread2("UKBB/coding609.tsv")
# df_self_reported <- fread2(csv, select = c(paste0("20002-0.", 0:33),
#                                            paste0("20002-1.", 0:33),
#                                            paste0("20002-2.", 0:33),
#                                            paste0("20002-3.", 0:33))) %>%
#   mutate_all(~ as.character(factor(., levels = coding609$coding,
#                                    labels = coding609$meaning)))

# Depression ever diagnosed by a professional
df_dep <- fread2(csv, select = c(paste0("20544-0.", 1:16))) %>%
  mutate_all(~ ifelse(. == 11, "F330", NA)) %>%
  select_if(is.character)

df_ICD9 <- fread2(csv, colClasses = "character", select = c(
  paste0("40013-", 0:14, ".0"),    # cancer
  paste0("41271-0.", 0:46)         # diagnosis
))

id_icd10_count <- bind_cols(df_ICD10, df_dep) %>%
  mutate_all(~ ifelse(. == "", NA, .)) %>%
  mutate(id = row_number()) %>%
  tidyr::pivot_longer(-id, values_to = "code", values_drop_na = TRUE) %>%
  group_by(id, code) %>%
  summarise(count = n(), .groups = "drop") %>%
  transmute(id, vocabulary_id = "ICD10", code, count)
str(id_icd10_count)

id_icd9_count <- df_ICD9 %>%
  mutate_all(~ ifelse(. == "", NA, .)) %>%
  mutate(id = row_number()) %>%
  tidyr::pivot_longer(-id, values_to = "code", values_drop_na = TRUE) %>%
  group_by(id, code) %>%
  summarise(count = n(), .groups = "drop") %>%
  transmute(id, vocabulary_id = "ICD9", code, count)
str(id_icd9_count)

# https://phewascatalog.org/phecodes
map_icd9 <- fread2("tmp-data/phecode_icd9_rolled.csv", colClasses = "character")
phecode_map_icd9 <- transmute(map_icd9, vocabulary_id = "ICD9",
                              code = ICD9, phecode = PheCode)


# https://github.com/PheWAS/PheWAS
library(PheWAS)
phecodes <- createPhenotypes(
  rbind(id_icd10_count, id_icd9_count),
  id.sex = data.frame(id = seq_along(sex), sex = c("F", "M")[sex + 1L]),
  vocabulary.map = mutate_at(rbind(phecode_map_icd10, phecode_map_icd9),
                             "code", ~ sub("\\.", "", .)),
  min.code.count = 1,
  add.phecode.exclusions = TRUE,
  full.population.ids = seq_along(sex)
)

saveRDS(phecodes[colSums(phecodes, na.rm = TRUE) > 500],
        "data/all_phecodes.rds")

# some verifs
table(pheno = phecodes$`185`, sex, exclude = NULL)
#        sex
# pheno        0      1   <NA>
#   FALSE      0 195543  13476
#   TRUE       0  10905    285
#   <NA>  264796  17019    481
table(pheno = phecodes$`174.1`, sex, exclude = NULL)
table(pheno = phecodes$`654.2`, sex, exclude = NULL)
table(pheno = phecodes$`296.2`, sex, exclude = NULL)
