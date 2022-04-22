#' The function converts file in VCF format from GATK into a data frame
#'
#' @param file The name of genotype file in VCF format from GATK SNP calling pipeline
#' @param BulkA The sample name of one bulk
#' @param BulkB The sample name of the other bulk
#' @param chromlist Default value is “NULL”.A string vector of the chromosomes to be used in the analysis, used to filter out unneeded contigs etc
#'
#' @return a data frame containing columns for Read depth (DP),Reference Allele Depth (AD_REF), Alternative Allele Depth (AD_ALT) and Genoytype Quality (GQ) for each bulk (indicated by .Bulk1 and .Bulk2 column name suffix)
#' @export vcf2table
#' @import dplyr
#' @import vcfR 
#' @importFrom readr parse_number
#' @importFrom tidyr %>% separate
#' 
 
vcf2table <- function(file,BulkA,BulkB,chromlist = NULL){
  cat("Let's Begin \n")
  vcf <- vcfR::read.vcfR(file)
  if(!is.null(vcf)){
    cat("import vcf file finished \n")
    vfix <- dplyr::as_tibble(vcf@fix[, c("CHROM", "POS", "REF", "ALT")]) %>% dplyr::mutate(Key = seq(1:nrow(.)))
    if(is.character(vfix$CHROM)==TRUE){
      vfix$CHROM <- parse_number(vfix$CHROM)
    }
    cat("The number of markers is",nrow(vfix),"\n")
    if((is.character(BulkA)==FALSE)|(is.character(BulkB)==FALSE)){
      message("Error! Our calculate STOP! Please check your parameter whether exactly right?")
    }else{
      vcount <- vcfR::extract_gt_tidy(vcf, format_fields = c("AD", "DP", "GQ"), gt_column_prepend = "", alleles = FALSE) %>%
        separate(
          col = "AD",
          into = c("AD.REF", "AD.ALT"),
          sep = ",",
          extra = "merge",
          convert = TRUE
        )
      cat("Extracting AD DP GQ finished \n")
      vcount1 <- vcount[vcount$Indiv==BulkA,]
      vcount1 <- select(vcount1,-Indiv) 
      vcount2 <- vcount[vcount$Indiv==BulkB,]
      vcount2 <- select(vcount2,-Indiv) 
      colnames(vcount1) <- c("Key",paste0(BulkA,".",names(vcount1[,2:ncol(vcount1)])))
      colnames(vcount2) <- c("Key",paste0(BulkB,".",names(vcount2[,2:ncol(vcount2)])))
      set <- dplyr::full_join(vfix,vcount1,by = "Key")
      data <- dplyr::full_join(set,vcount2,by = "Key") %>% dplyr::select(-Key)
      cat("Data integration completed \n")
      if (!is.null(chromlist)) {
        cat("Keep only wanted chromosomes",chromlist,"\n")
        data <- data[data$CHROM %in% chromlist, ]
      }
      cat("Finish! \n")
      return(data)
    }
  }else{
    message("Error! Our calculate STOP! Please check your vcf file connection exactly right?")
  }
}
