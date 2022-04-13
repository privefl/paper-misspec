# https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
pgen_zst <- runonce::download_file(
  "https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1",
  dir = "tmp-data")
pvar_zst <- runonce::download_file(
  "https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1",
  dir = "tmp-data")
psam <- runonce::download_file(
  "https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1",
  dir = "tmp-data", fname = "all_phase3.psam")

library(bigsnpr)
plink2 <- download_plink2("tmp-data")

sysglue <- function(...) system(glue::glue(...))
pgen <- sub('\\.zst', '', pgen_zst)
if (!file.exists(pgen)) sysglue("{plink2} --zst-decompress {pgen_zst} > {pgen}")
pvar <- sub('\\.zst', '', pvar_zst)
if (!file.exists(pvar)) sysglue("{plink2} --zst-decompress {pvar_zst} > {pvar}")

# https://dx.doi.org/10.1016%2Fj.ajhg.2017.03.004
outliers <- c("NA20314", "HG00731", "HG00732", "HG01880", "HG01882",
              "HG01944", "HG02497", "NA20320", "NA20321")
bigsnpr:::write.table2(cbind.data.frame(0, outliers), (tmp <- tempfile()))

sysglue("{plink2} --pfile tmp-data/all_phase3",
        " --remove {tmp} --keep-founders",
        " --max-alleles 2 --mac 10 --autosome",
        " --make-bed --out tmp-data/all_phase3",
        " --threads {nb_cores()[[1]]} --memory 50e3")

obj.bed <- bed("tmp-data/all_phase3.bed")

snp_info <- snp_attach("data/UKBB_HM3_val.rds")$map

library(dplyr)
matched <- snp_match(
  transmute(obj.bed$map, chr = as.integer(chromosome), pos = physical.pos,
            a0 = allele2, a1 = allele1, beta = 1),
  transmute(snp_info, chr = chromosome, pos = physical.pos, a0 = allele2, a1 = allele1)
) %>%
  arrange(`_NUM_ID_`)
# 23,675,459 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,054,315 variants have been matched; 0 were flipped and 1,054,315 were reversed.

all.equal(matched$`_NUM_ID_`, rows_along(snp_info))  # all matched

rel <- snp_plinkKINGQC(plink2, obj.bed$bedfile, thr.king = 2^-2.5,
                       make.bed = FALSE, ncores = nb_cores())

snp_readBed2(obj.bed$bedfile, backingfile = "data/1000G_HM3",
             ind.row = which(!obj.bed$fam$sample.ID %in% rel$IID2),
             ind.col = matched$`_NUM_ID_.ss`,
             ncores = nb_cores())

obj.bigsnp <- snp_attach("data/1000G_HM3.rds")
obj.bigsnp$fam <- dplyr::left_join(
  obj.bigsnp$fam,
  bigreadr::fread2("https://figshare.com/ndownloader/files/31080292"))
table(obj.bigsnp$fam$Population)
# ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD
#  91  53  86  93  99 103 105  94  99  99  91 101 113
# IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI
# 107 101 104  99  96  85  64  84  96 102  98 107 108
table(obj.bigsnp$fam$`Super Population`)
# AFR AMR EAS EUR SAS
# 645 344 504 503 482

snp_save(obj.bigsnp)
