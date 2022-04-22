#' The function to filter SNPs based on read depth and quality of each bulk
#'
#' @param file The data frame exported by vcf2table() or the data frame who has the same structure as the output of vcf2table()
#' @param Bulk.DP The minimum total read depth for a SNP
#' @param Bulk2.Ref The minimum read depth of a reference SNP allele 
#' @param Bulk2.Alt The minimum read depth of a alternative SNP allele
#' @param Depth.diff The maximum absolute differences between the read depths of the bulks
#' @param min.GQ The minimum Genotype Quality output by GATK
#' @param verbose Default value is “TRUE”. If the value is "TRUE", the number of SNPs filtered in each step will print
#'
#'
#' @return a data frame containing SNPs after filtering
#' @export BSA_filter
#'

                         
BSA_filter <- function(file,Bulk.DP,Bulk2.Ref,Bulk2.Alt,Depth.diff,min.GQ,verbose = TRUE){
  org_count <- nrow(file)
  now_count <- nrow(file)
  # Filter by sample DP
  if (!missing(Bulk.DP)) {
    if (verbose) {
      cat("Filtering by total sample read depth: Each sample DP >= ",Bulk.DP,"\n")
    }
    flag = which((file[,11] >= Bulk.DP) & (file[,7] >= Bulk.DP) == TRUE)
    file <- file[flag, ]
    if (verbose) {
      cat("Filtered ",now_count - nrow(file)," SNPs \n")
    }
    now_count <- nrow(file)
  }
  # Filter by sample Bulk2
  if (!missing(Bulk2.Ref)) {
    # Filter by Bulk2.Ref
    if (verbose) {
      cat("Filtering by Bulk2 referece allele read depth: Bulk2 AD.Ref >= ",Bulk2.Ref,"\n")
    }
    flag = which((file[,9]>= Bulk2.Ref) == TRUE)
    file <- file[flag, ]
    if (verbose) {
      cat("Filtered ",now_count - nrow(file)," SNPs \n")
    }
    now_count <- nrow(file)
  }
  # Filter by sample Bulk2
  if (!missing(Bulk2.Alt)) {
    # Filter by Bulk2.Alt
    if (verbose) {
      cat("Filtering by Bulk2 referece allele read depth: Bulk2 AD.Alt >= ",Bulk2.Alt,"\n")
    }
    flag = which((file[,10]>= Bulk2.Alt) == TRUE)
    file <- file[flag, ]
    if (verbose) {
      cat("Filtered ",now_count - nrow(file)," SNPs \n")
    }
    now_count <- nrow(file)
  }
  # Filter by difference between high and low bulks
  if (!missing(Depth.diff)) {
    if (verbose) {
      cat("Filtering by depth difference between bulks <= ",Depth.diff,"\n")
    }
      flag = which((abs(file[,7] - file[,11]) <= Depth.diff) == TRUE)
      file <- file[flag, ]
    if (verbose) {
      cat("Filtered ",now_count - nrow(file)," SNPs \n")
    }
    now_count <- nrow(file)
  }
  # Filter by Genotype Quality
  if (!missing(min.GQ)) {
      if (verbose) {
        cat("Filtering by Genotype Quality: GQ >= ",min.GQ,"\n")
      }
    flag = which(((file[,8]>= min.GQ) & (file[,12]>= min.GQ)) == TRUE)
    file <- file[flag, ]
      if (verbose) {
        cat("Filtered ",now_count - nrow(file)," SNPs \n")
      }
    now_count <- nrow(file)
      } 
  if (verbose) {
    cat("Original SNP number: ",org_count,", Filtered: ",org_count - now_count,", Remaining: ",now_count,"\n")
  }
  return(as.data.frame(file))
}
