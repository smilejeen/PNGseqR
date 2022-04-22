#' The function to perform BSA analysis with △SNP algorithm
#'
#' @param file The data frame exported by BSA_filter()  
#' @param Windowsize A numeric value indicates the window size (in base pairs) containing SNPs subjected to calculating △SNP
#' @param threshold From 0 to 1, default value is “0.95”. A numeric threshold value to extract candidate region by fractile quantile
#' @param pval A numeric threshold value, default value is “0.05”
#'
#'
#' @return a list containing three data frame,the data frame named "total.result" containing all the calculate result by △SNP algorithm,the data frame named "sigthreshold.result" containing the significant SNPs under fractile quantile set by parameter "threshold" calculated by △SNP algorithm,the data frame named "sigpval.result" containing the significant SNPs selected by P-value set by parameter "pval" calculated by G△SNP algorithmm
#' @export DeltaSNP
#' @importFrom stats predict
#' @importFrom locfit locfit lp
#' 
                         
DeltaSNP <- function(file,Windowsize,threshold,pval){
  cat("Let's Begin \n")
  if(!is.null(file)){
    if((is.numeric(Windowsize)==TRUE)&(Windowsize>=0)){
      file[,c(5:12)] <- as.numeric(unlist(file[,c(5:12)]))
      SNPindex.B1 <- file[,6]/file[,7]
      SNPindex.B2 <- file[,10]/file[,11]
      deltaSNP <- abs(SNPindex.B1 - SNPindex.B2)
      cat("Getting the delta SNP value \n")
      re <- cbind(file,SNPindex.B1,SNPindex.B2,deltaSNP)
      re <- re[!is.na(re$deltaSNP),]
      names(re) <- c(colnames(re)[1:12],"SI.Bulk1","SI.Bulk2","deltaSNP")
      cat("Getting the smooth result by Nadaraya Watson kernel regression\n")
      df <- data.frame()
      result <- data.frame()
      for (i in 1:length(unique(re$CHROM))) {
        dta <- re[re$CHROM==unique(re$CHROM)[i],]
        tricubeDS <- predict(locfit(deltaSNP~lp(POS,h=Windowsize,deg = 0),dta),dta,POS)
        dta <- cbind(dta,tricubeDS)
        result <- rbind(result,dta)
      }
      cat("Estimating the significance of each locus from smooth result by permutation test\n")
      mt <- mean(result$tricubeDS)
      if(length(result$tricubeDS)<=50000){
        sample <- sample(result$tricubeDS,999)
        sm <- sample - mt
        pvalue <- unlist(lapply(result$tricubeDS,function(x){(length(sm[sm>(x-mt)])+0.5)/1000}))
      }else{
        sample <- sample(result$tricubeDS,9999)
        sm <- sample - mt
        pvalue <- unlist(lapply(result$tricubeDS,function(x){(length(sm[sm>(x-mt)])+0.5)/10000}))
      }
      cat("Extracting the regions higher than threshold quantile\n")
      result <- cbind(result,pvalue)
      names(result) <- c(colnames(result)[1:16],"pvalue")
      threshold.value <- sort(as.numeric(result$tricubeDS))[1+(length(result$tricubeDS)-1)*threshold]
      sigthreshold.result <- result[result$tricubeDS>=threshold.value,]
      sigpval.result <- result[result$pvalue<=pval,]
      result <- list(result,sigthreshold.result,sigpval.result)
      names(result) <- c("total.result","sigthreshold.result","sigpval.result")
      cat("Finish! \n")
      return(result)
    }else{
      message("Error! Our calculate STOP! Please check your parameter whether exactly right?")
    }
  }else{
    message("Error! Our calculate STOP! Please check your import file whether exist?")
  }
}