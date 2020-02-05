## code to prepare the `nhrs` dataset is shown below

## nhrs is a table join between targets_and_families.csv,
## a file containing a list of all NHRs with suppl data,
## and biomart_nhr_ensg_hgnc.csv, a file created
## using biomart

library(tidyverse)
library(readxl)
library(usethis)
library(biomaRt)
library(janitor)
library(here)

if (!file.exists(here::here("data-raw/targets_and_families.csv"))) {
  download.file(url="http://www.guidetopharmacology.org/DATA/targets_and_families.csv",
                destfile=here::here("data-raw/targets_and_families.csv"))
}

raw <- read_csv(file=here::here("data-raw/targets_and_families.csv")) %>%
  janitor::clean_names()
nhr_db <- filter(raw, type == "nhr")


if (!file.exists(here::here("data-raw/biomart_nhr_ensg_hgnc.csv"))) {
  ensembl <- useMart(biomart="ensembl",
                     host="http://useast.ensembl.org",
                     dataset="hsapiens_gene_ensembl")
  bm <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                   mart=ensembl,
                   filters="entrezgene_id",
                   values=unique(nhr_db$human_entrez_gene) %>% as.character)
  write_csv(x=bm,
            path=here::here("data-raw", "biomart_nhr_ensg_hgnc.csv"),
            quote_escape=FALSE,
            col_names=TRUE
  )
}

biomart_nhr_ensg_hgnc <- read_csv(here::here("data-raw/biomart_nhr_ensg_hgnc.csv"))

nhr_db2 <- nhr_db %>%
  left_join(biomart_nhr_ensg_hgnc %>%
              distinct, by=c("human_entrez_gene"="entrezgene_id"))

nhrs <- nhr_db2 %>%
  dplyr::select("hgnc_symbol", "hgnc_id", "hgnc_name", "human_entrez_gene",
         "ensembl_gene_id", "synonyms") %>%
  rename(entrez_gene_id=human_entrez_gene)

write_csv(x=nhrs, path="data-raw/nhrs.csv")
usethis::use_data(nhrs, overwrite=TRUE, compress="xz")
