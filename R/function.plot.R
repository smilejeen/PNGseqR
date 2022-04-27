#' The function to output BSA analysis result into scatter and line plots
#'
#' @param file The data frame exported by the BSA analysis functions 
#' @param chromlist A string vector of the chromosome names used to plot line plots specific chromosomes, default value is "NULL"
#' @param type A character to indicate the selected algorithm
#' @param file.name A character to name the output plot
#' @param threshold From 0 to 1, default value is "NULL". The red threshold line displayed in the output plot
#' @param pval Default value is "NULL". A numeric threshold value drawn in blue line in the output plot

#'
#' @return a scatter plot for statistic value and a line plot for tricube-smoothed statistic value of whole genome, if parameter "chromlist" is not "NULL", will also output a series line plot for tricube-smoothed statistic value of each chromosome 
#' @export plot_BSA
#' @importFrom grDevices dev.off pdf png
#' @importFrom graphics abline axis par points title
#' 
 
plot_BSA <- function(file,chromlist=NULL,type=c("DeltaSNP","Bayesian","ED","Gprime"),file.name,threshold=NULL,pval=NULL){
  if(!is.null(file)){
      if(is.character(file.name)==TRUE){
        if((type=="DeltaSNP")|(type=="Bayesian")|(type=="ED")|(type=="Gprime")){
          cat("Let's Begin \n")
          chr.length <- c()
          chroms <- as.character(unique(data_filt$CHROM))
          file$POS <- as.numeric(file$POS)
          for (i in 1:length(chroms)) {
            dta <- file[file$CHROM==chroms[i],]
            chrl <- max(dta$POS)
            chr.length <- c(chr.length,chrl)
          }
          for(i.chr in chroms){
            if(i.chr==chroms[1]){
              pos.elment1 <- chr.length[1]
              pos.elment2 <- pos.elment1/2
              pos <- pos.elment1
              pos.chr <- pos.elment2
            }else{
              pos.elment2 <- pos.elment1+chr.length[which(chroms==i.chr)]/2
              pos.elment1 <- sum(chr.length[1:which(chroms==i.chr)])
              pos <- c(pos,pos.elment1)
              pos.chr <- c(pos.chr,pos.elment2)}
          }
          for(j in chroms[2:length(chroms)]){
            if(j==chroms[2]){
              chr.data <- file[file[,1]==j,]
              col2.data <- chr.data[,2]+pos[1]
              chr.data[,2] <-  col2.data
              mht.data <- rbind(file[file[,1]==chroms[1],],chr.data)
            }else{
              chr.data <- file[file[,1]==j,]
              col2.data <- chr.data[,2]+pos[which(chroms==j)-1]
              chr.data[,2] <- col2.data
              mht.data <- rbind(mht.data,chr.data)
            }
          }
          if(type == "Bayesian"){
            cat("Plotting Bayesian method Scatter diagram \n")
            png(paste0(file.name,".point.png"),height=480,width=1080) # 命名
            y <- mht.data[,15]
            x <- mht.data[,2]
            y.max <- ceiling(max(as.numeric(y),na.rm=TRUE))
            color.array <- rep(c("blue4","green4"), len = length(chroms))
            plot(x,y,type="p",cex=1,cex.axis=2,cex.lab=1.5,xlab="Chromosomes region",ylab="Probability of Linkage",xlim=c(0,pos[length(pos)]),
                 ylim=c(0,y.max),xaxs="i",yaxs="i",xaxt="n",font=2)
            for(k in chroms){
              x.zc <- mht.data[mht.data[,1]==k,2]
              y.zc <- mht.data[mht.data[,1]==k,15]
              points(x.zc,y.zc,col=color.array[which(chroms==k)],pch=20,cex=1)
            }
            axis(side = 1,at=pos.chr,labels=chroms,cex.axis=2,family="serif",font=2)
            par(family="serif",font.lab=2,font.axis=2)
            dev.off()
            cat("Bayesian method Scatter diagram finished \n")
            png(paste0(file.name,".smooth.png"),height=480,width=1080) # 命名
            y0 <- mht.data[,16]
            x0 <- mht.data[,2]
            y0.max <- max(y0,na.rm=TRUE)
            color.array <- rep(c("black"), len = length(chroms))
            cat("Plotting Bayesian method line Chart \n")
            plot(x0,y0,type="l",cex=1,cex.axis=2,cex.lab=1.5,xlab="Chromosomes region",ylab="Tricube smooth of Probability of Linkage",xlim=c(0,pos[length(pos)]),
                 ylim=c(0,y0.max),xaxs="i",yaxs="i",xaxt="n",font=2)
            if (!is.null(threshold)){
              threshold.value <- sort(as.numeric(file$tricubePL))[1+(length(file$tricubePL)-1)*threshold]
              sig.result <- file[file$tricubePL>=threshold.value,]
              h.v <- sig.result[which.min(sig.result$tricubePL),16]
              abline(h=h.v,lwd=3.5,col="red")
            }
            if (!is.null(pval)){
              p.result <- file[file$pvalue<=pval,]
              p.v <- p.result[which.min(p.result$tricubePL),16]
              abline(h=p.v,lwd=3.5,col="blue")
            }
            for(k in chroms){
              x0.zc <- mht.data[mht.data[,1]==k,2]
              y0.zc <- mht.data[mht.data[,1]==k,16]
              points(x0.zc,y0.zc,col=color.array[which(chroms==k)],pch=20,cex=1)
            }
            axis(side = 1,at=pos.chr,labels=chroms,cex.axis=2,family="serif",font=2)
            par(family="serif",font.lab=2,font.axis=2)
            dev.off()
            cat("Bayesian method line Chart finished \n")
            if (!is.null(chromlist)){
              cat("Plotting Bayesian method line Chart of selected chromosomes \n")
              for (i in chromlist) {
                png(paste0(file.name," Chromosome ",i,".smooth.png"),height=500,width=430) # 命名
                test <- file[file$CHROM==i,]
                x0 <- test$POS/1000000
                y0 <- test$tricubePL
                plot(x0,y0,type="l",cex=1.5,cex.axis=1.5,cex.lab=1.5,lwd = 4.5,xlab="Genomic Position",ylab="Tricube smooth of Bayesian",xlim=c(0,x0[length(x0)]),
                     ylim=c(0,max(y0)+1),xaxs="i",yaxs="i",font=2)
                title(main = paste0("chromosome ",i),cex.main=2)
                if (!is.null(threshold)){
                  threshold.value <- sort(as.numeric(file$tricubePL))[1+(length(file$tricubePL)-1)*threshold]
                  sig.result <- file[file$tricubePL>=threshold.value,]
                  h.v <- sig.result[which.min(sig.result$tricubePL),16]
                  abline(h=h.v,lwd=3.5,col="red")
                }
                if (!is.null(pval)){
                  p.result <- file[file$pvalue<=pval,]
                  p.v <- p.result[which.min(p.result$tricubePL),16]
                  abline(h=p.v,lwd=3.5,col="blue")
                }
                dev.off()
            }
              cat("Bayesian method line Chart of selected chromosomes finished \n")
            }
          }else if (type == "ED"){
            cat("Plotting ED method Scatter diagram \n")
            png(paste0(file.name,".point.png"),height=480,width=1080) # 命名
            y <- mht.data[,14]
            x <- mht.data[,2]
            y.max <- max(y,na.rm=TRUE)
            color.array <- rep(c("blue4","green4"), len = length(chroms))
            plot(x,y,type="p",cex=1,cex.axis=2,cex.lab=1.5,xlab="Chromosomes region",ylab="ED^Power",xlim=c(0,pos[length(pos)]),
                 ylim=c(0,y.max),xaxs="i",yaxs="i",xaxt="n",font=2)
            for(k in chroms){
              x.zc <- mht.data[mht.data[,1]==k,2]
              y.zc <- mht.data[mht.data[,1]==k,14]
              points(x.zc,y.zc,col=color.array[which(chroms==k)],pch=20,cex=1)
            }
            axis(side = 1,at=pos.chr,labels=chroms,cex.axis=2,family="serif",font=2)
            par(family="serif",font.lab=2,font.axis=2)
            dev.off()
            cat("ED method Scatter diagram finished \n")
            png(paste0(file.name,".smooth.png"),height=480,width=1080) # 命名
            y0 <- mht.data[,15]
            x0 <- mht.data[,2]
            y0.max <- max(y0,na.rm=TRUE)
            color.array <- rep(c("black"), len = length(chroms))
            cat("Plotting ED method line Chart \n")
            plot(x0,y0,type="l",cex=1,cex.axis=2,cex.lab=1.5,xlab="Chromosomes region",ylab="Tricube smooth of ED power",xlim=c(0,pos[length(pos)]),
                 ylim=c(0,y0.max),xaxs="i",yaxs="i",xaxt="n",font=2)
            if (!is.null(threshold)){
              threshold.value <- sort(as.numeric(file$tricubeED))[1+(length(file$tricubeED)-1)*threshold]
              sig.result <- file[file$tricubeED>=threshold.value,]
              h.v <- sig.result[which.min(sig.result$tricubeED),15]
              abline(h=h.v,lwd=3.5,col="red")
            }
            if (!is.null(pval)){
              p.result <- file[file$pvalue<=pval,]
              p.v <- p.result[which.min(p.result$tricubeED),15]
              abline(h=p.v,lwd=3.5,col="blue")
            }
            for(k in chroms){
              x0.zc <- mht.data[mht.data[,1]==k,2]
              y0.zc <- mht.data[mht.data[,1]==k,15]
              points(x0.zc,y0.zc,col=color.array[which(chroms==k)],pch=20,cex=1)
            }
            axis(side = 1,at=pos.chr,labels=chroms,cex.axis=2,family="serif",font=2)
            par(family="serif",font.lab=2,font.axis=2)
            dev.off()
            cat("ED method line Chart finished \n")
            if (!is.null(chromlist)){
            cat("Plotting ED method line Chart of selected chromosomes \n")
            for (i in chromlist) {
              png(paste0(file.name," Chromosome ",i,".smooth.png"),height=500,width=430) # 命名
              test <- file[file$CHROM==i,]
              x0 <- test$POS/1000000
              y0 <- test$tricubeED
              plot(x0,y0,type="l",cex=1.5,cex.axis=1.5,cex.lab=1.5,lwd = 4.5,xlab="Genomic Position",ylab="Tricube smooth of ED",xlim=c(0,x0[length(x0)]),
                   ylim=c(0,max(y0)+1),xaxs="i",yaxs="i",font=2)
              title(main = paste0("chromosome ",i),cex.main=2)
              if (!is.null(threshold)){
                threshold.value <- sort(as.numeric(file$tricubeED))[1+(length(file$tricubeED)-1)*threshold]
                sig.result <- file[file$tricubeED>=threshold.value,]
                h.v <- sig.result[which.min(sig.result$tricubeED),14]
                abline(h=h.v,lwd=3.5,col="red")
              }
              if (!is.null(pval)){
                p.result <- file[file$pvalue<=pval,]
                p.v <- p.result[which.min(p.result$tricubeED),14]
                abline(h=p.v,lwd=3.5,col="blue")
              }
              dev.off()
            }
              cat("ED method line Chart of selected chromosomes finished \n")
            }
          } else if (type == "Gprime"){
            cat("Plotting Gprime method Scatter diagram \n")
            png(paste0(file.name,".point.png"),height=480,width=1080) # 命名
            y <- mht.data[,13]
            x <- mht.data[,2]
            y.max <- max(y,na.rm=TRUE)
            color.array <- rep(c("blue4","green4"), len = length(chroms))
            plot(x,y,type="p",cex=1,cex.axis=2,cex.lab=1.5,xlab="Chromosomes region",ylab="G-value",xlim=c(0,pos[length(pos)]),
                 ylim=c(0,y.max),xaxs="i",yaxs="i",xaxt="n",font=2)
            for(k in chroms){
              x.zc <- mht.data[mht.data[,1]==k,2]
              y.zc <- mht.data[mht.data[,1]==k,13]
              points(x.zc,y.zc,col=color.array[which(chroms==k)],pch=20,cex=1)
            }
            axis(side = 1,at=pos.chr,labels=chroms,cex.axis=2,family="serif",font=2)
            par(family="serif",font.lab=2,font.axis=2)
            dev.off()
            cat("Gprime method Scatter diagram finished \n")
            png(paste0(file.name,".smooth.png"),height=480,width=1080) # 命名
            y0 <- mht.data[,14]
            x0 <- mht.data[,2]
            y0.max <- max(y0,na.rm=TRUE)
            color.array <- rep(c("black"), len = length(chroms))
            cat("Plotting Gprime method line Chart \n")
            plot(x0,y0,type="l",cex=1,cex.axis=2,cex.lab=1.5,xlab="Chromosomes region",ylab="Tricube smooth of G-value",xlim=c(0,pos[length(pos)]),
                 ylim=c(0,y0.max),xaxs="i",yaxs="i",xaxt="n",font=2)
            if (!is.null(threshold)){
              threshold.value <- sort(as.numeric(file$tricubeG))[1+(length(file$tricubeG)-1)*threshold]
              sig.result <- file[file$tricubeG>=threshold.value,]
              h.v <- sig.result[which.min(sig.result$tricubeG),14]
              abline(h=h.v,lwd=3.5,col="red")
            }
            if (!is.null(pval)){
              p.result <- file[file$pvalue<=pval,]
              p.v <- p.result[which.min(p.result$tricubeG),14]
              abline(h=p.v,lwd=3.5,col="blue")
            }
            for(k in chroms){
              x0.zc <- mht.data[mht.data[,1]==k,2]
              y0.zc <- mht.data[mht.data[,1]==k,14]
              points(x0.zc,y0.zc,col=color.array[which(chroms==k)],pch=20,cex=1)
            }
            axis(side = 1,at=pos.chr,labels=chroms,cex.axis=2,family="serif",font=2)
            par(family="serif",font.lab=2,font.axis=2)
            dev.off()
            cat("Gprime method line Chart finished \n")
            if (!is.null(chromlist)){
            cat("Plotting Gprime method line Chart of selected chromosomes \n")
            for (i in chromlist) {
              png(paste0(file.name," Chromosome ",i,".smooth.png"),height=500,width=430) # 命名
              test <- file[file$CHROM==i,]
              x0 <- test$POS/1000000
              y0 <- test$tricubeG
              plot(x0,y0,type="l",cex=1.5,cex.axis=1.5,cex.lab=1.5,lwd = 4.5,xlab="Genomic Position",ylab="Tricube smooth of G-value",xlim=c(0,x0[length(x0)]),
                   ylim=c(0,max(y0)+1),xaxs="i",yaxs="i",font=2)
              title(main = paste0("chromosome ",i),cex.main=2)
              if (!is.null(threshold)){
                threshold.value <- sort(as.numeric(file$tricubeG))[1+(length(file$tricubeG)-1)*threshold]
                sig.result <- file[file$tricubeG>=threshold.value,]
                h.v <- sig.result[which.min(sig.result$tricubeG),14]
                abline(h=h.v,lwd=3.5,col="red")
              }
              if (!is.null(pval)){
                p.result <- file[file$pvalue<=pval,]
                p.v <- p.result[which.min(p.result$tricubeG),14]
                abline(h=p.v,lwd=3.5,col="blue")
              }
              dev.off()
            }
              cat("Gprime method line Chart of selected chromosomes finished \n")
            }
          } else if (type == "DeltaSNP"){
            cat("Plotting deltaSNP method Scatter diagram \n")
            png(paste0(file.name,".point.png"),height=480,width=1080) # 命名
            y <- mht.data[,15]
            x <- mht.data[,2]
            y.max <- max(y,na.rm=TRUE)
            color.array <- rep(c("blue4","green4"), len = length(chroms))
            plot(x,y,type="p",cex=1,cex.axis=2,cex.lab=1.5,xlab="Chromosomes region",ylab="deltaSNP",xlim=c(0,pos[length(pos)]),
                 ylim=c(0,y.max),xaxs="i",yaxs="i",xaxt="n",font=2)
            for(k in chroms){
              x.zc <- mht.data[mht.data[,1]==k,2]
              y.zc <- mht.data[mht.data[,1]==k,15]
              points(x.zc,y.zc,col=color.array[which(chroms==k)],pch=20,cex=1)
            }
            axis(side = 1,at=pos.chr,labels=chroms,cex.axis=2,family="serif",font=2)
            par(family="serif",font.lab=2,font.axis=2)
            dev.off()
            cat("deltaSNP method Scatter diagram finished \n")
            png(paste0(file.name,".smooth.png"),height=480,width=1080) # 命名
            y0 <- mht.data[,16]
            x0 <- mht.data[,2]
            y0.max <- max(y0,na.rm=TRUE)
            color.array <- rep(c("black"), len = length(chroms))
            cat("Plotting deltaSNP method line Chart \n")
            plot(x0,y0,type="l",cex=1,cex.axis=2,cex.lab=1.5,xlab="Chromosomes region",ylab="Tricube smooth of deltaSNP",xlim=c(0,pos[length(pos)]),
                 ylim=c(0,y0.max),xaxs="i",yaxs="i",xaxt="n",font=2)
            if (!is.null(threshold)){
              threshold.value <- sort(as.numeric(file$tricubeDS))[1+(length(file$tricubeDS)-1)*threshold]
              sig.result <- file[file$tricubeDS>=threshold.value,]
              h.v <- sig.result[which.min(sig.result$tricubeDS),16]
              abline(h=h.v,lwd=3.5,col="red")
            }
            if (!is.null(pval)){
              p.result <- file[file$pvalue<=pval,]
              p.v <- p.result[which.min(p.result$tricubeDS),16]
              abline(h=p.v,lwd=3.5,col="blue")
            } 
            for(k in chroms){
              x0.zc <- mht.data[mht.data[,1]==k,2]
              y0.zc <- mht.data[mht.data[,1]==k,16]
              points(x0.zc,y0.zc,col=color.array[which(chroms==k)],pch=20,cex=1)
            }
            axis(side = 1,at=pos.chr,labels=chroms,cex.axis=2,family="serif",font=2)
            par(family="serif",font.lab=2,font.axis=2)
            dev.off()
            cat("DeltaSNP method line Chart finished \n")
            if (!is.null(chromlist)){
            cat("Plotting DeltaSNP method line Chart of selected chromosomes \n")
            for (i in chromlist) {
              png(paste0(file.name," Chromosome ",i,".smooth.png"),height=500,width=430) # 命名
              test <- file[file$CHROM==i,]
              x0 <- test$POS/1000000
              y0 <- test$tricubeDS
              plot(x0,y0,type="l",cex=1.5,cex.axis=1.5,cex.lab=1.5,lwd = 4.5,xlab="Genomic Position",ylab="Tricube smooth of deltaSNP",xlim=c(0,x0[length(x0)]),
                   ylim=c(0,max(y0)+1),xaxs="i",yaxs="i",font=2)
              title(main = paste0("chromosome ",i),cex.main=2)
              if (!is.null(threshold)){
                threshold.value <- sort(as.numeric(file$tricubeDS))[1+(length(file$tricubeDS)-1)*threshold]
                sig.result <- file[file$tricubeDS>=threshold.value,]
                h.v <- sig.result[which.min(sig.result$tricubeDS),14]
                abline(h=h.v,lwd=3.5,col="red")
              }
              if (!is.null(pval)){
                p.result <- file[file$pvalue<=pval,]
                p.v <- p.result[which.min(p.result$tricubeDS),14]
                abline(h=p.v,lwd=3.5,col="blue")
              }
              dev.off()
            }
              cat("DeltaSNP method line Chart of selected chromosomes finished \n")
            }
          }
        }else{
          message("Error! Our calculate STOP! Please check your parameter type whether exactly right?")
        }
      }else{
        message("Error! Our calculate STOP! Please check your parameter file.name whether exactly right?")
      }
  }else{
    message("Error! Our calculate STOP! Please check your import file whether exist?")
  }
}

