#' Table of Nuclear Hormone Receptors (NHRs)
#'
#' This object contains a table of all known NHRs and was adapted from the
#' 'guidetopharmacology' website (see references). It was joined with a
#' bioMart table to include ensemble gene ids, which are commonly used
#' gene symbols.
#' @format A tibble with 54 rows and 6 variables:
#' \describe{
#'   \item{hgnc_symbol}{the HUGO gene nomenclature committee (HGNC) symbol
#'   (letters and numbers, ex. THRB)}
#'   \item{hgnc_id}{the HUGO gene nomenclature committee (HGNC) symbol
#'   (a number, ex. 11799)}
#'   \item{hgnc_name}{the HUGO gene nomenclature committee (HGNC) gene name
#'   (ex. "Thyroid hormone receptor beta")}
#'   \item{entrez_gene_id}{the entrez gene id (a number, ex. 7068)}
#'   \item{ensembl_gene_id}{the ensembl gene id (ex. ENSG00000151090, always
#'   starts with ENSG)}
#'   \item{synonyms}{words or gene symbols in the literature that refer to
#'   the same gene}
#' }
#' @source \url{
#' http://www.guidetopharmacology.org/DATA/targets_and_families.csv}
#' @source \url{http://www.biomart.org/}
"nhrs"
