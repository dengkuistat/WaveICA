\name{correlation}
\alias{correlation}
\title{Computing correlations between any two samples.}
\usage{
 correlation(data_before,data_after,method="pearson")

)
}
\arguments{
\item{data_before}{Sample-by-matrix original metabolomics data.}
\item{data_after}{Sample-by-matrix  metabolomics data after applying batch removal methods.}
\item{method}{Specify the method to compute correlation, one of the pearson and spearman.}

}
\value{
A matrix containing the correlation coefficients between any two samples.
}

\description{
Computing correlations between any two samples.
}
\author{
Kui deng
\email{dengkui_stat@163.com}
}


