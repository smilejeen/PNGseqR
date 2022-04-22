#' The function to extract genes fall in the candidate region and judge whether gene is differential expression
#' 
#' @param BSA_data A data frame contains three columns, including the chromosome name, the start and the end positions of the candidate regions
#' @param gtf_data The name of reference genome in GTF format
#' @param DEG_data A vector containing the names of differential expression genes,default is "NULL"
#'
#' @return if parameter "DEG_data" is set "NULL" will return a data frame contain the genes fall in the candidate region,if paramater "DEG_data" is a vector provide the gene id of differential expression will return a list, the first data frame in the list contains the genes fall in the candidate region,the second data frame contains the differential expression genes fall in the candidate region
#' @export result_anno
#' @importFrom rtracklayer import
#' @importFrom readr parse_number
#' @importFrom dplyr filter left_join
#'
                                    
result_anno <- function(BSA_data,
                       gtf_data,
                       DEG_data=NULL){
  cat("Let's begin \n")
  cat("Import gtf file \n")
  gtf <- import(gtf_data,"gtf")
  cat("Import gtf finished \n")
  gtf <- as.data.frame(gtf[gtf$type == "gene",])
  cat("So there are",nrow(gtf),"genes \n" )
  colsremain <- c("start","end","gene_id")
  df <- gtf[,colnames(gtf)%in%colsremain]
  if(is.character(gtf[,1])==TRUE){
    new_gtf <- cbind(parse_number(gtf[,1]),df)
  } else {
    new_gtf <- cbind(gtf[,1],df)
  }
  names(new_gtf)<-c("CHROM","start","end","gene_id")
  names(BSA_data) <- c("CHROM","start","end")
  cat("Exracting genes in candidate region \n")
  re <- data.frame()
  for (i in 1:dim(BSA_data)[1]) {
    cf <- BSA_data[i,1]
    sf <- BSA_data[i,2]
    ef <- BSA_data[i,3]
    df <- filter(new_gtf,CHROM==cf,start<=ef,end>=sf)
    re <- rbind(re,df)
  }
   result <- re
   if (!is.null(DEG_data)) {
     cat("Selecting DEG genes in candidate region \n")
     DEG_id <- DEG_data$id
     sig_result <- re[re$gene_id %in% DEG_id,]
     deg_result <- left_join(sig_result,DEG_data,by=c('gene_id'='id'))
     result <- list(re,deg_result)
     names(result) <- c("genes_in_region","DEGgenes_in_region")
   }
  cat("Finished! \n")
   return(result)
}



