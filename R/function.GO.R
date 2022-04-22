#' The function to perform GO analysis
#'
#' @param GO_term A data frame contains at least two columns, one column is genes' name and the other column is the GO term
#' @param DEG_id A vector containing the names of the differential expression genes
#' @param type  A character indicates the mode of GO analysis. There are four options: "CC" indicates cellular component; "BP" indicates biological process, "MF" indicates molecular function, and "all" indicates all three modes. The default is “all” 
#' @param algorithm A character indicates the algorithm for GO analysis. Three options "classic","weight01", and "elim" indicate three different algorithms
#' @param statistic A character indicates two different significance test methods. Two options "fisher" and "ks" indicates two different methods
#' @param top_num A numeric value to indicate the number of significant GO terms
#' @param plot Default value is “TRUE”. If "TRUE", the acyclic, bubble, and histogram plots will be created
#'
#'
#' @return a data frame of the gene ontology analysis result and plot directed acyclic graphs, bubble graphs, histograms if parameter plot is set "TRUE"
#' @export GO_analysis
#' @importFrom utils write.table
#' @importFrom methods new
#' @importFrom grDevices dev.off pdf png
#' @import topGO
#' @import ggplot2
#' 
                               
GO_analysis <- function(GO_term,
                        DEG_id,
                        type=c("BP","MF","CC","all"),
                        algorithm=c("classic","weight01","elim"),
                        statistic=c("fisher","ks"),
                        top_num,
                        plot = TRUE){
  cat("Let's begin \n")
  cat("So the number of differential expression gene is ",length(DEG_id),"\n")
  names(go_term) <- c("id","GO_terms")
  go_split <- split(go_term$GO_terms, go_term$id)
  cat("So the number of Species gene is ",length(go_split),"\n")
  re <- NULL
  id <- names(go_split)
  cat("Building gene2go mapping,please waiting!","\n")
  for (i in 1:length(go_split)) {
    term <- paste(go_split[[i]],collapse=",")
    re <- rbind(re,term)
  }
  cat("Finish bulding gene2go mapping!","\n")
  result <- cbind(id,re)
  write.table(result,"gene2go.map",sep = "\t",quote = F,row.names = F,col.names = F)
  gene2go <- topGO::readMappings(file = "gene2go.map")
  cat("Prepare GO matrix","\n")
  genenames <- names(gene2go)
  genelist <- factor(as.integer(genenames %in% DEG_id)) 
  # 这里会生成一个factor，有两个levels：0和1，其中1表示感兴趣的基因。
  names(genelist) <- genenames
  if(type!="all"){
  GOdata <- new("topGOdata", ontology=type, allGenes = genelist,annot = annFUN.gene2GO, gene2GO = gene2go)
  cat("GO analysing!","\n")
  node <-  topGO::runTest(GOdata, algorithm = algorithm, statistic = statistic)
  if(statistic=="ks"){
    table <- topGO::GenTable(GOdata,KS = node,ranksOf = algorithm,topNodes = attributes(node)$geneData[4])
  }else{
    table <- topGO::GenTable(GOdata,Fis = node,ranksOf = algorithm,topNodes = attributes(node)$geneData[4])
  }
  Gene.ratio <- table$Significant/table$Annotated
  Category <- rep(i,dim(table)[1])
  result <- cbind(table,Gene.ratio,Category)
  result <- result[which((result$Gene.ratio)>0),]
  result <- result[1:top_num,]
  names(result) <- c("GO.ID","Term","Annotated","Significant","Expected","Pvalue","Gene.ratio","Category")
  result$Pvalue <- as.numeric(result$Pvalue)
  if (plot){
  cat("Plotting Directed acyclic graph \n")
  pdf(paste0(type,"_DAG.pdf"),width = 8,height = 6.5)
  topGO::showSigOfNodes(GOdata, score(node), firstSigNodes = 5, useInfo = "all")
  dev.off()
  cat("Plotting Bubble graph \n")
  pdf(paste0(type,"_bubble.pdf"),height=6,width=8)
  print(ggplot2::ggplot(result,aes(Gene.ratio,Term))+ggplot2::geom_point(aes(size = Significant,color = Pvalue))+
          ggplot2::scale_color_gradient(low = "blue",high = "red")+ggplot2::theme_bw()+
          ggplot2::labs(color=expression(-log[10](Pvalue)),size="Count",x="Gene Ratio",y="GO Terms"))
  dev.off()
  cat("Plotting histogram \n")
  pdf(paste0(type,"_histogram.pdf"),height=6,width=8)
  print(ggplot2::ggplot(result,aes(y=Significant,x=Term,fill=Category)) + ggplot2::geom_bar(stat="identity",position = "dodge") + ggplot2::facet_grid(Category~.,scales = "free",space = "free") + ggplot2::coord_flip() + ggplot2::theme_bw()+ ggplot2::labs(x="Terms",y="Significant Count"))
  dev.off()
  }
  cat("Finish! \n")
  return(result)
  }else{
    type <- c("BP","MF","CC")
    result <- data.frame()
    for (i in type) {
      GOdata <- new("topGOdata", ontology=i, allGenes = genelist,annot = annFUN.gene2GO, gene2GO = gene2go)
      cat("GO analysing for ",i," !","\n")
      node <-  topGO::runTest(GOdata, algorithm = algorithm, statistic = statistic)
      if(statistic=="ks"){
        table <- topGO::GenTable(GOdata,KS = node,ranksOf = algorithm,topNodes = attributes(node)$geneData[4])
      }else{
        table <- topGO::GenTable(GOdata,Fis = node,ranksOf = algorithm,topNodes = attributes(node)$geneData[4])
      }
      Gene.ratio <- table$Significant/table$Annotated
      Category <- rep(i,dim(table)[1])
      table <- cbind(table,Gene.ratio,Category)
      table <- table[which((table$Gene.ratio)>0),]
      table <- table[1:top_num,]
      result <- rbind(result,table)
      cat(i,"analysis finish! \n")
      if (plot){
        cat("Plotting",i," Directed acyclic graph \n")
        pdf(paste0(i,"_DAG.pdf"),width = 8,height = 6.5)
        topGO::showSigOfNodes(GOdata, score(node), firstSigNodes = 5, useInfo = "all")
        dev.off()
      }
    }
    names(result) <- c("GO.ID","Term","Annotated","Significant","Expected","Pvalue","Gene.ratio","Category")
    result$Pvalue <- as.numeric(result$Pvalue)
    if (plot){
      cat("Plotting Bubble graph \n")
      pdf("all_bubble.pdf",height=6,width=8)
      print(ggplot2::ggplot(result,aes(Gene.ratio,Term))+ggplot2::geom_point(aes(size = Significant,color = Pvalue))+
              ggplot2::scale_color_gradient(low = "blue",high = "red")+ggplot2::theme_bw()+
              ggplot2::labs(color=expression(-log[10](Pvalue)),size="Count",x="Gene Ratio",y="GO Terms"))
      dev.off()
      cat("Plotting histogram \n")
      pdf("all_histogram.pdf",height=6,width=8)
      print(ggplot2::ggplot(result,aes(y=Significant,x=Term,fill=Category)) + ggplot2::geom_bar(stat="identity",position = "dodge") + ggplot2::facet_grid(Category~.,scales = "free",space = "free") + ggplot2::coord_flip() + ggplot2::theme_bw() + ggplot2::labs(x="Terms",y="Significant Count"))
      dev.off()
    }
    cat("Finish! \n")
    return(result)
  }
}
