#' @title correlation
#' @description Computing correlations between any two samples.
#' @author Kui Deng
#' \email{dengkui_stat@163.com}
#' @param data_before Sample-by-matrix original metabolomics data.
#' @param data_after Sample-by-matrix  metabolomics data after applying batch removal methods.
#' @param method Specify the method to compute correlation, one of the pearson and spearman.
#' @return A matrix containing the correlation coefficients between any two samples.
#' @export

correlation<-function(data_before,data_after,method="pearson"){
  data_before_cor<-cor(t(data_before),method=method)
  data_after_cor<-cor(t(data_after),method=method)

    C<-matrix(NA,(ncol(data_before_cor)*ncol(data_before_cor)-ncol(data_before_cor))/2,4)
  count<-0
  for(i in 1:(ncol(data_before_cor)-1)){
    for(j in (i+1):ncol(data_before_cor)){
      count=count+1
      C[count,3]<-data_before_cor[i,j]
      C[count,4]<-data_after_cor[i,j]
      C[count,1]<-i
      C[count,2]<-j
    }
  }

  C<-as.data.frame(C)
  colnames(C)<-c("var_before","var_after","cor_before","cor_after")
  return(C)
}

