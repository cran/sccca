#' Performs fisher exact test to get the significant overlap between genes for cell type assignment
#'
#' This Function is used to perform fisher exact test to get cell types
#'
#' @usage fisher_test(ref, gene_overlap)
#'
#' @param ref reference gene set.
#'
#' @param gene_overlap genes that pass the correlation threshold.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return vector of p-value and overlap.
#'
#' @export
#'
#' @examples fisher_test(c("PAX8","PAX6","TP53","AOC3","LIPF"), c("LIPF","PAX8","PAX6","TP53","TSHB","AOC3"))
#'
fisher_test <- function(ref, gene_overlap){
  # Calculate the size of each set
  n_pd <- length(ref)
  n_ad <- length(gene_overlap)

  # Calculate the size of the overlap (intersection) between the sets
  overlap <- length(intersect(ref, gene_overlap))

  # Perform Fisher's Exact Test
  fisher_test_result <- fisher.test(matrix(c(overlap, n_pd - overlap, n_ad - overlap, 0), nrow = 2),
                                     simulate.p.value = T)
  return(c(fisher_test_result$p.value, overlap))
}
#' Performs phyper test to get the significant overlap between genes for cell type assignment
#'
#' This Function is used to perform phyper test to get cell types
#'
#' @usage phyper_test(ref, overlap)
#'
#' @param ref reference gene set.
#'
#' @param overlap genes that pass the correlation threshold.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return vector of p-value and overlap.
#'
#' @export
#'
#' @examples phyper_test(c("PAX8","PAX6","TP53","AOC3","LIPF"), c("LIPF","PAX8","PAX6","TP53","TSHB","AOC3"))
#'
phyper_test <- function(ref, overlap) {
  # Convert vectors to sets to get unique elements
  set1 <- unique(ref)
  set2 <- unique(overlap)
  # Calculate the size of the intersection
  intersection_size <- length(intersect(set1, set2))
  # Calculate the total number of unique elements
  total_elements <- length(union(set1, set2))
  # Total number of elements in the first vector
  total_elements1 <- length(set1)
  # Total number of elements in the second vector
  total_elements2 <- length(set2)
  # Calculate the p-value using phyper
  p_value <- phyper(intersection_size - 1, total_elements1, total_elements1, total_elements2, lower.tail = FALSE)
  # Return the p-value
  return(c(p_value, intersection_size))
}
#' Performs parallel function on two lists
#'
#' This Function is used to perform parallel function on two lists
#'
#' @usage enrich_genes(ref_list, overlap_list, func)
#'
#' @param ref_list reference list.
#'
#' @param overlap_list overlap list.
#'
#' @param func function to be applied.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#' 
#' @return list where each element is the result of applying the function `func` to the corresponding elements of `ref_list` and `overlap_list`. 
#'
#' @export
#'
enrich_genes <- function(ref_list, overlap_list, func) {
  mapply(func, ref_list, overlap_list, SIMPLIFY = FALSE)
}
