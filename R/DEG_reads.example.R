#' Reads count example data
#' 
#' A dataset containing two repeat reads count for mutant pool and wild type pool (calculated from NCBI, accession number: PRJNA699154) 
#'
#' @docType data
#' 
#' @usage data(DEG_reads.example)
#' 
#' 
#' @format Dataframe, one row per gene, with five columns
#' \describe{
#' \item{Gene id}{Gene symbol/name used for differential expression analysis refers to. Data type: character}
#' \item{MU1}{Reads number of each gene in the first mutant repeat. Data type: integer}
#' \item{MU2}{Reads number of each gene in the second mutant repeat. Data type: integer}
#' \item{WT1}{Reads number of each gene in the first wild type repeat. Data type: integer}
#' \item{WT2}{Reads number of each gene in the second wild type repeat. Data type: integer}
#' }
#' 
"DEG_reads.example"