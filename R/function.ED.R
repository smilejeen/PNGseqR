#' The function to perform BSA analysis with Euclidean distance algorithm
#'
#' @param file The data frame exported by BSA_filter()  
#' @param power A numeric value indicates the power to raise ED. It was suggested to be 4~6,default value is "4"
#' @param Windowsize A numeric value indicates the window size (in base pairs) containing SNPs subjected to calculating Euclidean distance
#' @param threshold From 0 to 1, default value is "0.95". A numeric threshold value to extract candidate region by fractile quantile
#' @param pval A numeric threshold value, default value is "0.05"
#'
#'
#' @return a list containing three data frame,the data frame named "total.result" containing all the calculate result by Euclidean Distance algorithm,the data frame named "sigthreshold.result" containing the significant SNPs under fractile quantile set by parameter "threshold" calculated by Euclidean Distance algorithm,the data frame named "sigpval.result" containing the significant SNPs selected by P-value set by parameter "pval" calculated by Euclidean Distance algorithm 
#' @export ED
#' @importFrom stats predict
#' @importFrom locfit locfit lp
#' 

ED <- function(file,power = 4,Windowsize,threshold = 0.95,pval = 0.05){
  cat("Let's Begin \n")
  if(!is.null(file)){
    if((is.numeric(power)==TRUE)&(Windowsize>=0)){
      file[,c(5:12)] <- as.numeric(unlist(file[,c(5:12)]))
      B1.REF.REQ = file[,5]/file[,7]
      B2.REF.REQ = file[,9]/file[,11]
      B1.ALT.REQ = file[,6]/file[,7]
      B2.ALT.REQ = file[,10]/file[,11]
      ED <- sqrt((B1.REF.REQ-B2.REF.REQ)^2+(B1.ALT.REQ-B2.ALT.REQ)^2)
      ED.power <- ED^power
      cat("Getting the calculate result by Euclidean distance analysis\n")
      data.new <- cbind(file,ED,ED.power)
      re <- data.new[!is.na(data.new$ED.power),]
      names(re) <- c(colnames(re)[1:12],"ED","ED.power")
      cat("Getting the smooth result by Nadaraya Watson kernel regression\n")
      df <- data.frame()
      result <- data.frame()
      for (i in 1:length(unique(re$CHROM))) {
        dta <- re[re$CHROM==unique(re$CHROM)[i],]
        tricubeED <- predict(locfit(ED.power~lp(POS,h=Windowsize,deg = 0),dta),dta,POS)
        dta <- cbind(dta,tricubeED)
        result <- rbind(result,dta)
      }
      cat("Estimating the significance of each locus from smooth result by permutation test\n")
      mt <- mean(result$tricubeED)
      if(length(result$tricubeED)<=50000){
        sample <- sample(result$tricubeED,999)
        sm <- sample - mt
        pvalue <- unlist(lapply(result$tricubeED,function(x){(length(sm[sm>(x-mt)])+0.5)/1000}))
      }else{
        sample <- sample(result$tricubeED,9999)
        sm <- sample - mt
        pvalue <- unlist(lapply(result$tricubeED,function(x){(length(sm[sm>(x-mt)])+0.5)/10000}))
      }
      cat("Extracting the regions higher than threshold quantile\n")
      result <- cbind(result,pvalue)
      names(result) <- c(colnames(result)[1:15],"pvalue")
      threshold.value <- sort(as.numeric(result$tricubeED))[1+(length(result$tricubeED)-1)*threshold]
      sigthreshold.result <- result[result$tricubeED>=threshold.value,]
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
