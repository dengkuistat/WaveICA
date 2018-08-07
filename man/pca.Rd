\name{pca}
\alias{pca}
\title{PCA score plot.}
\usage{
 pca(
  data,
  id,
  plot,
  col=c("red","blue","black","green"),
  pch=c(19,17,15,18),
  label=c("Control","Disease","QC")
)
}
\arguments{
\item{data}{Sample-by-matrix metabolomics data.}
\item{id}{group.}
\item{plot}{1: two dimension PCA score plots; 2: three dimension PCA score plot.}
\item{col}{color of plot}
\item{pch}{shape of plot.}
\item{label}{label of plot}

}
\value{
PCA score plot.
}

\description{
Plot PCA score plot.
}
\author{
Kui deng
\email{dengkui_stat@163.com}
}


