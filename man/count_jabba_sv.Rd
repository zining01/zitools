% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{count_jabba_sv}
\alias{count_jabba_sv}
\title{count_jabba_sv}
\usage{
count_jabba_sv(jabba_rds = NA_character_, id = "", verbose = FALSE)
}
\arguments{
\item{jabba_rds}{(character)}

\item{id}{(character) default ""}

\item{verbose}{(logical) default FALSE}
}
\value{
data.table with columns
- pair
- n.jun
- n.loose
}
\description{
Given a JaBbA graph, count how many loose ends and junctions
}
