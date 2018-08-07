\name{SVM_MODEL}
\alias{SVM_MODEL}
\title{k-fold cross validation of SVM}
\usage{
SVM_MODEL(X,
          Y,
          kernel="radial",
          kfold=5,
          seed=1234)
}
\arguments{
\item{data}{Sample-by-matrix metabolomics data.}
\item{Y}{group of the samples.}
\item{kernel}{kernel function, the default is radial.}
\item{kfold}{Specify the fold of the cross validation, the default is 5-fold cross validation.}
\item{seed}{seed to repeat the process.}

}
\value{
A list containing the AUC value of predictive values.}

\description{
Cross validation of SVM model.}
\author{
Kui deng
\email{dengkui_stat@163.com}
}


