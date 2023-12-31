% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{cncp_analysis}
\alias{cncp_analysis}
\title{cncp_analysis}
\usage{
cncp_analysis(
  seg = NULL,
  gs = "~/projects/JaBbA_pipeline/db/hcc1954.sv.gs.rds",
  method = "jabba",
  id = "sample",
  max.dist = 10000,
  mask = "~/projects/gGnome/files/zc_stash/maskA_re.rds",
  field = "cn",
  simplify = FALSE,
  loose.only = FALSE,
  ignore.strand = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{seg}{(character) path to seg file (can be BigWig, GRanges, or gGraph)}

\item{gs}{(character) path to gold standard SVs (as GRangesList)}

\item{method}{(character) analysis method (one of jabba, ascat, sequenza, default jabba)}

\item{id}{(character) sample ID}

\item{max.dist}{(numeric) bp threshold under which a called CNCP is considered to overlap a gold-standard CNCP}

\item{mask}{(character) masked ranges (ignore breakpoints in this region)}

\item{simplify}{(logical) simplify the input segments when calling CNCP?}

\item{loose.only}{(logical) only consider the loose ends}

\item{verbose}{(logical) default FALSE}
}
\value{
list with entries:
- $cncp.gr: GRanges of CN change points in the segmentation + distance to nearest gold-standard point
- $cncp.res: containing F1 score, precision, and recall of CN change points
}
\description{
This is a wrapper function for CN change point analysis.
It compares a CN segmentation of a genome to a set of gold-standard break points.
}
