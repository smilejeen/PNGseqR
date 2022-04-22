#' Genotype data
#' 
#' A dataset containing SNP information in  Variant Call Format(VCF) 
#'
#' @docType data
#' 
#' @usage data(vcf_data.example)
#' 
#' 
#' @format vcf, one row per SNP, with eleven columns
#' \describe{
#' \item{CHROM}{The chromosome which SNP fall in. Data type: character}
#' \item{POS}{The position on each chromosome of SNP, in base pairs. Data type: character}
#' \item{ID}{The name of each SNP. Data type: character}
#' \item{REF}{The base type at this position of reference genome. Data type: character}
#' \item{ALT}{Compare to the reference genome,the of alternative base type at this position. Data type: character}
#' \item{QUAL}{Quality value of variation site detection. Data type: character}
#' \item{FILTER}{Used to judge if the site passed the filter criteria. Data type: character}
#' \item{FORMAT}{Including the genotype information(GT),allele depth in each sample information(AD), total depth of this site information(DP), Genotype quality value(GQ), Likelihood value of genotype(PL). Data type: character}
#' \item{B73_MU}{The value cprresponding to FORMAT in B73_MU sample. Data type: character}
#' \item{B73_WT}{The value cprresponding to FORMAT in B73_WT sample. Data type: character}
#' }
#' 
"vcf_data.example"