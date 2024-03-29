% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat.R
\name{enrich_genes}
\alias{enrich_genes}
\title{Performs parallel function on two lists}
\usage{
enrich_genes(ref_list, overlap_list, func)
}
\arguments{
\item{ref_list}{reference list.}

\item{overlap_list}{overlap list.}

\item{func}{function to be applied.}
}
\value{
list where each element is the result of applying the function `func` to the corresponding elements of `ref_list` and `overlap_list`.
}
\description{
This Function is used to perform parallel function on two lists
}
\author{
Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
}
