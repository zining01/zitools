% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{reads2clouds}
\alias{reads2clouds}
\title{reads2clouds}
\usage{
reads2clouds(reads, thresh = 5000, read.size = 101)
}
\arguments{
\item{reads}{(data.table) must have columns seqnames, start, BX}

\item{thresh}{(numeric) distance between barcodes for merging to clouds}

\item{read.size}{(numeric) base pairs of average read?}
}
\value{
data.table representing read clouds
modified version of reads to clouds
}
\description{
reads2clouds
}
