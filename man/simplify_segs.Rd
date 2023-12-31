% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{simplify_segs}
\alias{simplify_segs}
\title{simplify_segs}
\usage{
simplify_segs(gr = NULL, sort.gr = TRUE, verbose = FALSE)
}
\arguments{
\item{gr}{(GRanges)}

\item{sort.gr}{(logical) default TRUE}

\item{verbose}{(logical) default FALSE}
}
\value{
simplifed GRanges with NA segments removed
}
\description{
simplify CN segmentation of GRanges object. Must have field $cn.
}
