% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calculate_cor_mat}
\alias{calculate_cor_mat}
\title{Performs aggregation based on cell clusters and condition. Then, it calculates correlation matrix of genes}
\usage{
calculate_cor_mat(expression_mat, condition  = NULL, clusters, assay = "RNA")
}
\arguments{
\item{expression_mat}{Seurat object that contains the expression matrix.}

\item{condition}{column name of the condition in th meta data of the Seurat object.}

\item{clusters}{column name of the cluster numbers in the meta data of the Seurat object.}

\item{assay}{the assay to be used default is set to RNA}
}
\value{
correlation matrix of genes.
}
\description{
This Function is used to perform cell aggregation by averaging the expression of scRNA-seq matrix and then perform correlation matrix
}
\author{
Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
}
