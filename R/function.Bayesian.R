#' The function to perform BSA analysis with empirical Bayesian algorithm 
#'
#' @param file The data frame exported by BSA_filter()
#' @param Bulksize The number of samples in each pool
#' @param Genomelength The total length of reference genome
#' @param Windowsize A numeric value indicates the window size (in base pairs) containing SNPs subjected to calculating empirical Bayesian
#' @param threshold A numeric threshold value. Default value is "0.95"
#' @param pval A numeric threshold value. Default value is "0.05"
#'
#'
#' @return a list containing three data frame,the data frame named "total.result" containing all the calculate result by empirical Bayesian algorithm,the data frame named "sigthreshold.result" containing the significant SNPs under fractile quantile set by parameter "threshold" calculated by empirical Bayesian algorithm,the data frame named "sigpval.result" containing the significant SNPs selected by P-value set by parameter "pval" calculated by empirical Bayesian algorithm 
#' @export Bayesian
#' @importFrom stats dbinom predict
#' @importFrom locfit locfit lp
#'

                      
Bayesian <- function(file,Bulksize,Genomelength,Windowsize,threshold = 0.95,pval = 0.05){
  # calculate function
  calFUN <- function(x) {
    mu.ad=sum(x[1:2])
    wt.ad=sum(x[3:4])
    srof=which.min(x[1:2])
    ms=(x[1:2])[srof]
    ws=(x[3:4])[srof]
    if(ms==0){
      p1=mean(dbinom(ws,wt.ad,h.theta[h.theta>=0.5]))
      p2=mean(dbinom(ws,wt.ad,h.theta))
      pmt=(1/(1+mean((1-h.theta)^mu.ad+h.theta^mu.ad)*((1-p.theta)/p.theta)))
      pwt=p1/p2/2
    }
    else{
      pmt=0
      pwt=0
    }
    ppp=pmt*pwt
    c(pmt,pwt,ppp)
  }
  # filter data
  cat("Let's Begin \n")
  if(!is.null(file)){
    if((Windowsize>=0)&(Bulksize>=0)&(Genomelength>=0)){
      file[,c(5:12)] <- as.numeric(unlist(file[,c(5:12)]))
      d <- c(seq(0,Genomelength/100,by=0.01))/100 #为什么是0.01
      qd <- (1-.5*(1-exp(-2*d)))^(2*Bulksize)
      p.theta=sum((0.01/Genomelength)*(qd[-length(qd)]+qd[-1]))
      h.theta=c(file[,9]/apply(file[,9:10],1,sum),1-(file[,9]/apply(file[,9:10],1,sum)))
      dta <- as.matrix(file[,c(5,6,9,10)])
      cat("Getting the probabilities result by bayesian analysis\n")
      re <- cbind(file,t(apply(dta, 1, calFUN)))
      names(re) <- c(colnames(re)[1:12],"PMU","PWT","PPP")
      cat("Getting the smooth result by Nadaraya Watson kernel regression\n")
      df <- data.frame()
      result <- data.frame()
      for (i in 1:length(unique(re$CHROM))) {
        dta <- re[re$CHROM==unique(re$CHROM)[i],]
        tricubePL <- predict(locfit(PPP~lp(POS,h=Windowsize,deg = 0),dta),dta,POS)
        dta <- cbind(dta,tricubePL)
        result <- rbind(result,dta)
      }
      cat("Estimating the significance of each locus from smooth result by permutation test\n")
      mt <- mean(result$tricubePL)
      if(length(result$tricubePL)<=50000){
        sample <- sample(result$tricubePL,999)
        sm <- sample - mt
        pvalue <- unlist(lapply(result$tricubePL,function(x){(length(sm[sm>(x-mt)])+0.5)/1000}))
      }else{
        sample <- sample(result$tricubePL,9999)
        sm <- sample - mt
        pvalue <- unlist(lapply(result$tricubePL,function(x){(length(sm[sm>(x-mt)])+0.5)/10000}))
      }
      cat("Extracting the regions higher than threshold quantile\n")
      result <- cbind(result,pvalue)
      names(result) <- c(colnames(result)[1:16],"pvalue")
      threshold.value <- sort(as.numeric(result$tricubePL))[1+(length(result$tricubePL)-1)*threshold]
      sigthreshold.result <- result[result$tricubePL>=threshold.value,]
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