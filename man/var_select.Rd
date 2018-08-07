\name{var_select}
\alias{var_select}
\title{Variable selection}
\usage{
 var_select(data,label,t.test=T,Wilcox=T,AUC=T,FDR=T,VIP=T,FC=T,comps=3)

}
\arguments{
\item{data}{Sample-by-matrix metabolomics data.}
\item{label}{group of the samples.}
\item{t.test}{Whether to conduct t test.}
\item{Wilcox}{Whether to conduct wilcoxon test.}
\item{AUC}{Whether to calculate AUC value.}
\item{FDR}{Whether to calculate FDR value.}
\item{VIP}{Whether to calculate VIP value through PLS analysis.}
\item{FC}{Whether to calculate fold change value.}
\item{comps}{When doing PLS analysis, specify the components, the default is 3.}

}
\value{
Data matrix of the results of statistical analysis that you specify.
}

\description{
Variable selection by different statistical analysis.
}
\author{
Kui deng
\email{dengkui_stat@163.com}
}


