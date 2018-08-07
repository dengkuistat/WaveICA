#' @title SVM_MODEL
#' @description Cross validation of SVM model.
#' @author Kui Deng
#' \email{dengkui_stat@163.com}
#' @param data Sample-by-matrix metabolomics data.
#' @param Y group of samples.
#' @param kernel kernel function, the default is radial.
#' @param kfold Specify the fold of the cross validation, the default is 5-fold cross validation.
#' @param seed seed.
#' @return A list containing the AUC value of predictive values.
#' @export

SVM_MODEL<-function(X,
                    Y,
                    kernel="radial",
                    kfold=5,
                    seed=1234){
  set.seed(seed)
  kfold10 <- sample(rep(1:kfold,ceiling(nrow(X)/kfold))[1:nrow(X)])
  prob <- matrix(NA,nrow(X),1)
  for (i in (1:kfold)){
    cat(paste("##########fold",i,"##############\n"))
    library(e1071)
    library(pROC)
    svmmodel<-svm(x=as.matrix(X[kfold10!=i,,drop=F]),y=factor(Y[kfold10!=i]),kernel=kernel,probability = TRUE)
    prob[kfold10==i,1]<-attr(predict(svmmodel,as.matrix(X[kfold10==i,,drop=F]),probability=T),"probabilities")[,2]


  }
  roc<-roc(Y,as.numeric(prob))
  auc<-roc$auc
  result<-list(prob=prob,auc=auc)
}
