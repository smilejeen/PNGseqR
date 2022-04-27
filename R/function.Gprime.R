#' The function to perform BSA analysis with G-test algorithm
#'
#' @param file The data frame exported by BSA_filter()  
#' @param Windowsize A numeric value indicates the window size (in base pairs) containing SNPs subjected to calculating G statistic
#' @param threshold From 0 to 1, default value is “0.95”. A numeric threshold value to extract candidate region by fractile quantile
#' @param pval A numeric threshold value, default value is “0.05”
#'
#'
#' @return a list containing three data frame,the data frame named "total.result" containing all the calculate result by Gtest algorithm,the data frame named "sigthreshold.result" containing the significant SNPs under fractile quantile set by parameter "threshold" calculated by Gtest algorithm,the data frame named "sigpval.result" containing the significant SNPs selected by P-value set by parameter "pval" calculated by Gtest algorithm 
#' @export Gprime
#' @importFrom stats predict
#' @importFrom locfit locfit lp
#' 

Gprime <- function(file,Windowsize,threshold = 0.95,pval = 0.05){
  cat("Let's Begin \n")
  if(!is.null(file)){
    if((is.numeric(Windowsize)==TRUE)&(Windowsize>=0)){
      file[,c(5:12)] <- as.numeric(unlist(file[,c(5:12)]))
      exp <- cbind((file[,5]+file[,9])*(file[,11])/(file[,7]+file[,11]),(file[,5]+file[,9])*(file[,7])/(file[,7]+file[,11]),(file[,11])*(file[,6]+file[,10])/(file[,7]+file[,11]),(file[,6]+file[,10])*(file[,5]+file[,6])/(file[,7]+file[,11]))
      obs <- cbind(file[,9],file[,5],file[,10],file[,6])
      cat("Calculating G values from reads count matrix\n")
      G <-2*(apply((obs*log(obs/exp)),1,sum))
      data <- cbind(file,G)
      cat("Getting the smooth result by Nadaraya Watson kernel regression\n")
      result <- data.frame()
      for (i in 1:length(unique(file$CHROM))) {
        dta <- data[data$CHROM==unique(data$CHROM)[i],]
        tricubeG <- predict(locfit(G~lp(POS,h=Windowsize,deg = 0),dta),dta,POS)
        dta <- cbind(dta,tricubeG)
        result <- rbind(result,dta)
      }
      cat("Estimating the significance of each locus from smooth result by permutation test\n")
      mt <- mean(result$tricubeG)
      if(length(result$tricubeG)<=50000){
        sample <- sample(result$tricubeG,999)
        sm <- sample - mt
        pvalue <- unlist(lapply(result$tricubeG,function(x){(length(sm[sm>(x-mt)])+0.5)/1000}))
      }else{
        sample <- sample(result$tricubeG,9999)
        sm <- sample - mt
        pvalue <- unlist(lapply(result$tricubeG,function(x){(length(sm[sm>(x-mt)])+0.5)/10000}))
      }
      cat("Extracting the regions higher than threshold quantile\n")
      result <- cbind(result,pvalue)
      names(result) <- c(colnames(result)[1:14],"pvalue")
      threshold.value <- sort(as.numeric(result$tricubeG))[1+(length(result$tricubeG)-1)*threshold]
      sigthreshold.result <- result[result$tricubeG>=threshold.value,]
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