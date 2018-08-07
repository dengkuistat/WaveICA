#' @title var_select
#' @description Variable selection by different statistical analysis.
#' @author Kui Deng
#' \email{dengkui_stat@163.com}
#' @param data Sample-by-matrix metabolomics data.
#' @param label group of the samples.
#' @param t.test Whether to conduct t test.
#' @param Wilcox Whether to conduct wilcoxon test.
#' @param AUC Whether to calculate AUC value.
#' @param FDR Whether to calculate FDR value.
#' @param VIP Whether to calculate VIP value through PLS analysis.
#' @param FC Whether to calculate fold change value.
#' @param comps When doing PLS analysis, specify the components, the default is 3.
#' @return Data matrix of the results of statistical analysis that you specify.
#' @export

var_select<-function(data,label,t.test=T,Wilcox=T,AUC=T,FDR=T,VIP=T,FC=T,comps=3){
  library(pROC)
  library(parallel)
  data.result<-matrix(NA,ncol(data),6)
  label<-as.factor(label)
  data<-as.data.frame(data)
  if (length(unique(label))>2){
    stop("label must have two levels")
  }

  if (t.test==T){
    cat("########### ttest ###########\n")
    data.result[,1]<-as.numeric(mclapply(data,function(x) t.test(x~label)$p.value))
  }

  if (Wilcox==T){
    cat("########### wilcox ###########\n")
    data.result[,2]<-as.numeric(mclapply(data,function(x) wilcox.test(x~label)$p.value))
  }

  if (AUC==T){
    cat("########### AUC ###########\n")
    data.result[,3]<-as.numeric(mclapply(data,function(x) roc(label,x)$auc))
  }

  if (VIP==T){
    cat("########### VIP ###########\n")
    label1=c()
    for(i in 1:length(label)){
      if(label[i]==unique(label)[1]){label1[i]=unique(label)[2]}
      if(label[i]==unique(label)[2]){label1[i]=unique(label)[1]}
    }
    label12=cbind(label,label1)
    library(plsdepot)
    result<-plsreg2(data,label12,comps =comps,crosval=TRUE)
    data.result[,4]<-result$VIP[,comps]
  }

  if (FDR==T){
    library(fdrtool)
    cat("########### FDR ###########\n")
    p.value<-as.numeric(mclapply(data,function(x) t.test(x~label)$p.value))
    data.result[,5][!is.na(p.value)]<-fdrtool(p.value,statistic="pvalue",plot=F)$qval

  }
  if (FC==T){
    cat("########### Fold Change ###########\n")
    data.result[,6]<-apply(data,2,function(x) log(mean(x[label==unique(label)[2]])/mean(x[label==unique(label)[1]])))

  }
  data.result<-as.data.frame(data.result)
  rownames(data.result)<-colnames(data)
  colnames(data.result)<-c("P_t.test","P_wilcox","AUC","VIP","P.FDR","FC")
  return(data.result)
}

