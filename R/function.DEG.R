#' The function to perform differential expression analysis
#'
#' @param file A data frame containing read count of each pool
#' @param rep.num A numeric value indicates the number of replications in differential expression analysis
#' @param sig.level  numeric value indicates the FDR threshold in differential expression analysis 
#' @param file.name A character to name the output plot
#' @param exp.fold A numeric value indicates the log2FoldChange threshold in differential expression analysis 
#' @param plot Default value is "TRUE". If "TRUE", the volcano plot will be created. 
#'
#'
#' @return a data frame of the differential expression analysis result and a volcano plot of analysis result if parameter "plot" is set "TRUE" 
#' @export DEG_analysis
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom grDevices dev.off pdf png
#' @import ggplot2
#' 
                   
DEG_analysis <- function(file,rep.num,sig.level,file.name,exp.fold,plot = TRUE){
  if((is.data.frame(file)==TRUE)&(is.numeric(rep.num)==TRUE)&(is.numeric(sig.level)==TRUE)&(is.numeric(exp.fold)==TRUE)){
    cat("Let's Begin \n")
    count <- file
    cat("Import readscount files \n")
    levels <- c("B1","B2")
    names(count) <- c("id",paste0(rep(levels[1],rep.num),"-",1:rep.num),paste0(rep(levels[2],rep.num),"-",1:rep.num))
    rownames(count) <- count[,1]
    id <- count$id[1:(dim(count)[1])]
    count <- count[1:(dim(count)[1]),-1]
    readscount <- cbind(id,count)
    cat("Build matrix for differetial expression analysis \n")
    condition <- factor(c(rep(levels[1],rep.num),rep(levels[2],rep.num)),levels)
    colData <- data.frame(row.names=colnames(count),condition)
    dds <- DESeqDataSetFromMatrix(count, colData, design= ~ condition)
    cat("Perform differetial expression analysis \n")
    dds <- DESeq(dds)
    expdiff <- results(dds,contrast=c("condition",levels[1],levels[2]))
    expdiff <- cbind(id,expdiff)
    expdiff <- as.data.frame(expdiff)
    expdiff <- subset(expdiff,log2FoldChange!="NA") 
    expdiff <- subset(expdiff,padj!="NA") 
    sig.expdiff <- expdiff[((expdiff$padj <= sig.level)&(abs(expdiff$log2FoldChange) >= exp.fold)),]
    if (plot){
      cat("Plotting differential expression volcano diagram \n")
      png(paste0(file.name,"_volcano_plot.png"),height=540,width=480) # 命名
      threshold <- as.factor(ifelse(expdiff$padj <= sig.level & abs(expdiff$log2FoldChange) >= exp.fold ,ifelse(expdiff$log2FoldChange >= exp.fold ,'Up','Down'),'Not'))
      print(ggplot2::ggplot(expdiff,aes(x=log2FoldChange,y=-log10(padj),colour=threshold))+ggplot2::xlab("log2(Fold Change)")+ggplot2::ylab("-log10(qvalue)") + ggplot2::geom_point(size = 2,alpha=1)+ggplot2::ylim(0,20)+ggplot2::xlim(min(expdiff$log2FoldChange,na.rm=TRUE),max(expdiff$log2FoldChange,na.rm=TRUE))+ggplot2::scale_color_manual(values=c("blue","grey", "red"))+ggplot2::geom_vline(xintercept = c(-exp.fold, exp.fold), lty = 2,colour="#000000")+ggplot2::geom_hline(yintercept = c(-log10(sig.level)), lty = 2,colour="#000000")+ggplot2::theme(axis.text=element_text(size=20),axis.title=element_text(size=20))) 
      dev.off()
    }
    result <- list(expdiff,sig.expdiff)
    names(result) <- c("DEG_total","DEG_sig")
    return(result)
  }else{
    message("Error! Our calculate STOP! Please check your parameter whether exactly right?")
  }
}
  

