#' @title PCA
#' @description Plot PCA score plot.
#' @author Kui Deng
#' \email{dengkui_stat@163.com}
#' @param data Sample-by-matrix metabolomics data.
#' @param id group.
#' @param plot 1= 2 dimension PCA score plots; 2=3 dimension PCA score plot.
#' @param col color of plot.
#' @param pch shape of plot.
#' @param label label of plot.
#' @return PCA score plot.
#' @export

pca <- function(
  data,  ## Sample-by-variable matrix
  id,  ## group
  plot,  ## 1=2 dimension or 2=3 dimension PCA score plot
  col=c("red","blue","black","green"),  ## Color
  pch=c(19,17,15,18),  ## shape
  label=c("Control","Disease","QC")

)
{
  pc1 <- prcomp(data,retx=T,scale.=T)
  proportion<-pc1$sdev^2/sum(pc1$sdev^2)*100
  x<-pc1$x
  group <- id

  if(plot==1){
    plot(x[,1],x[,2],xlab=paste("PC1","(",round(proportion[1],digits =1),"%",")",sep=""),ylab=paste("PC2","(",round(proportion[2],digits =1),"%",")",sep=""),pch=pch[as.numeric(group)],col=col[as.integer(group)])
    #text(pc1$x[,1:2],row.names(data))
    legend("topright", legend=label,col=col,pch=pch, lwd=2)
  }



  if(plot==2){
    library(scatterplot3d)
    scatterplot3d(x[,1],x[,2],x[,3],angle=30,xlab=paste("PC1","(",round(proportion[1],digits =1),"%",")",sep=""),ylab=paste("PC2","(",round(proportion[2],digits =1),"%",")",sep=""),zlab=paste("PC3","(",round(proportion[3],digits =1),"%",")",sep=""),pch=pch[as.numeric(group)],color=col[as.integer(group)],box=FALSE)
    legend("topright", legend=label,col=col,pch=pch, lwd=2)
  }
}



