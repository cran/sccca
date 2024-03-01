#' Process the cell markers that overlap between the cell markers and scRNA matrix
#'
#' This Function is used to return the cell markers that overlap between the cell markers and scRNA matrix
#'
#' @usage filter_list(gene_list, passed_cells)
#'
#' @param gene_list list of unique genes of cell types.
#'
#' @param passed_cells cells types that pass the specified threshold.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return list of cell types which genes are found in the input matrix.
#'
#' @export
#'
filter_list <- function(gene_list, passed_cells) {
  filtered <- lapply(passed_cells, function(vec) vec[vec %in% gene_list])
  return(filtered[lengths(filtered) > 0])
}
#' Process the cell markers database and return the processed list
#'
#' This Function is used to process the cell markers database and return the processed list
#'
#' @usage process_markers(markers_df)
#'
#' @param markers_df data frame with markers named as gene_original and cell names as cell type.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return list of lists of the processed markers
#'
#' @export
#'
process_markers <- function(markers_df){

  gspositive_lst <- toupper(markers_df$gene_original)
  names(gspositive_lst) <- markers_df$celltype
  # filter by removing missing values
  gspositive_lst <- gspositive_lst[which(!is.na(gspositive_lst))]
  # merge named vectors with duplicate names
  gspositive_combined <- tapply(unlist(gspositive_lst, use.names = FALSE),
                                rep(names(gspositive_lst), lengths(gspositive_lst)), FUN = c)
  gspositive_combined  <- sapply(gspositive_combined,function(x) unique(x))
  return(gspositive_combined)
}

#' Process the cell markers that pass specific threshold in the gene correlation matrix
#'
#' This Function is used to return the cell markers that pass specific threshold in the gene correlation matrix
#'
#' @usage match_characters(genes, gene_mat)
#'
#' @param genes list of unique genes of cell types.
#'
#' @param gene_mat correlation matrix of genes.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return vector of genes names which overlap with the correlation matrix.
#'
#' @export
#'
match_characters <- function(genes, gene_mat) {
  match_genes <- genes %in% rownames(gene_mat)
  return(match_genes)
}
#' Calculate cell scores based on number of genes
#'
#' This Function is used to calculate cell scores based on number of genes
#'
#' @usage calculate_normalized_ratio(vec)
#'
#' @param vec list of genes of cell types.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return vector of cell scores based on the number of overlapped genes with the input matrix.
#'
#' @export
#'
calculate_normalized_ratio <- function(vec) {
  sum(vec) / length(vec)
}
#' Performs aggregation based on cell clusters and condition. Then, it calculates correlation matrix of genes
#'
#' This Function is used to perform cell aggregation by averaging the expression of scRNA-seq matrix and then perform correlation matrix
#'
#' @usage calculate_cor_mat(expression_mat, condition  = NULL, clusters, assay = "RNA")
#'
#' @param expression_mat Seurat object that contains the expression matrix.
#'
#' @param condition column name of the condition in th meta data of the Seurat object.
#'
#' @param clusters column name of the cluster numbers in the meta data of the Seurat object.
#'
#' @param assay the assay to be used default is set to RNA
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return correlation matrix of genes.
#'
#' @export
#'
calculate_cor_mat <- function(expression_mat, condition  = NULL, clusters, assay = "RNA"){
  # Get the average expression
  message("Perform Pesudobulk Aggregaion...")
  if (!is.null(condition)){
    expression_mat_avg <- AverageExpression(object = expression_mat, group.by = c(condition, clusters), assays = assay)
  }else{
    expression_mat_avg <- AverageExpression(object = expression_mat, group.by =  clusters, assays = assay)
  }
  message("Calculate Correlation matrix...")
  expression_mat_avg_mat <- as.matrix(expression_mat_avg[[1]])
  # Calculate the expression matrix
  cor_mat <- abs(cor(t(expression_mat_avg_mat)))
  return(cor_mat)
}
#' Filter the genes based on specific correlation threshold
#'
#' This Function is used to filter the gene correlation matrix based on user-defined threshold
#'
#' @usage filter_correlation(cor_mat, gene_list, threshold = 0.7)
#'
#' @param cor_mat correlation matrix generated from calculate_cor_mat function.
#'
#' @param gene_list cell markers that passed threshold.
#'
#' @param threshold absolute correlation threshold.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return vector of gene names that pass user-defined correlation threshold.
#'
#' @export
#'
filter_correlation <- function(cor_mat, gene_list, threshold = 0.7){
  gene_cor_mat <- cor_mat[rownames(cor_mat) %in% gene_list, colnames(cor_mat) %in% gene_list]
  diag(gene_cor_mat) <- NA

  genes_quantile <- apply(gene_cor_mat, 2, quantile, na.rm = T)
  passed_genes <- colnames(genes_quantile)[genes_quantile[3,] > threshold]
  return(passed_genes)
}
