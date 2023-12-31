% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{compute_cn_correlation}
\alias{compute_cn_correlation}
\title{compute_cn_correlation}
\usage{
compute_cn_correlation(
  gr = NULL,
  method = "cn",
  gs = "gs",
  min.width = 1000,
  verbose = FALSE
)
}
\arguments{
\item{gr}{(GRanges)}

\item{method}{(character) metadata field of method results}

\item{gs}{(character) metadata field of gold standard results}

\item{min.width}{(numeric) ignore segs shorter than this}

\item{verbose}{(logical)}
}
\value{
data.table with columns pearson.cn and spearman.cn
}
\description{
compute Pearson and Spearman correlation and RMSE
}
