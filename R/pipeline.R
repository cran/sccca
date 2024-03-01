#' Gets the associated cell types using correlation-based approach
#'
#' This Function is used to get the associated cell clusters using correlation-based approach
#'
#' @usage process_clus(cluster,sobj,assay="RNA",clus,markers,cor_m,m_t=0.9,c_t=0.7,test="p")
#'
#' @param cluster associated cluster name.
#'
#' @param sobj Seurat object.
#'
#' @param assay assay to be used default is set to RNA.
#'
#' @param clus cell clusters.
#'
#' @param markers cell markers database.
#'
#' @param cor_m gene correlation matrix.
#'
#' @param m_t overlap threshold between cell markers and expression matrix.
#'
#' @param c_t correlation threshold between genes.
#'
#' @param test statistical test that check if overlap is significant could be "p" for phyper or "f" for fisher.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return data frame of proposed cell types.
#'
#' @export
#'
process_clus <- function(cluster, sobj, assay = "RNA", clus, markers, cor_m, m_t = 0.9, c_t = 0.7, test = "p") {

  cell_clusters <- clus[clus[,1] == cluster, , drop = FALSE]
  # Get the matrix data associated with the cluster
  cluster_data <- sobj[,rownames(cell_clusters)]
  cluster_mat <- as.matrix(cluster_data@assays[[assay]]$counts)
  rownames(cluster_mat) <- toupper(rownames(cluster_mat))
  cluster_mat <- cluster_mat[apply(cluster_mat, 1, function(x) !all(x == 0)),]
  # Get the overlap between the matrix and the markers
  overlap_flag <- lapply(markers, match_characters, cluster_mat)
  # Get the normalized cell ratios
  normalized_cell_ratios <- data.frame(score = t(data.frame(lapply(overlap_flag, calculate_normalized_ratio), check.names = FALSE)))
  # Get cell types greater than threshold
  normalized_cell_ratios <- normalized_cell_ratios[normalized_cell_ratios$score >= m_t, , drop = FALSE]
  # Get cell types that have specific threshold
  passed_overlap <- markers[rownames(normalized_cell_ratios)]
  # Filter the markers based on correlation
  passed_overlap <- lapply(passed_overlap, function(cell_type_genes) {
    cell_type_mat <- filter_correlation(cor_m, cell_type_genes, c_t)
    return(cell_type_mat)
  })
  passed_overlap <- Filter(length, passed_overlap)
  # Get original list based on filtered markers
  ref_overlap <- markers[names(passed_overlap)]
  # Count gene occurrences and filter out unique genes
  genes_count <- table(unlist(ref_overlap))
  genes_unique <- names(genes_count[genes_count == 1])
  # Filter original and passed overlap lists based on unique genes
  ref_overlap <- filter_list(genes_unique, ref_overlap)
  passed_overlap <- filter_list(genes_unique, passed_overlap)
  ref_overlap <- ref_overlap[names(passed_overlap)]

  if(test == "f"){
    cell_types <- data.frame(t(data.frame(enrich_genes(ref_overlap, passed_overlap, fisher_test)))) %>%
      top_n(3) %>% mutate(cluster = cluster)
  }else {
    cell_types <- data.frame(t(data.frame(enrich_genes(ref_overlap, passed_overlap, phyper_test)))) %>%
      top_n(3) %>% mutate(cluster = cluster)
  }


  cell_types$cell_type <- rownames(cell_types)
  rownames(cell_types) <- NULL
  colnames(cell_types) <- c("p.value", "intersection_size", "cluster", "cell_type")
  cell_types <- cell_types[order(cell_types$cluster),]
  return(cell_types)
}

#' Run the pipeline for the cell type assignment
#'
#' This Function is used to run the main pipeline that does the cell type assignment
#'
#' @usage sccca(sobj,assay="RNA",cluster,marker,tissue,tt="a",cond,m_t=0.9,c_t=0.7,test="p",org="a")
#'
#' @param sobj Seurat object.
#'
#' @param assay assay to be used default is set to RNA.
#'
#' @param cluster colname in the mata.data that have the cell cluster numbers.
#'
#' @param marker cell markers database path.
#'
#' @param tissue specified tissue from which the data comes.
#' 
#' @param tt tissue type whether 'a' for all types 'n' for normal tissues only or "c" for cancer tissues.
#'
#' @param cond colname in the meta.data that have the condition names.
#'
#' @param m_t overlap threshold between cell markers and expression matrix.
#'
#' @param c_t correlation threshold between genes.
#'
#' @param test statistical test that check if overlap is significant could be "p" for phyper or "f" for fisher.
#'
#' @param org organism to be used that can be 'h' for human, 'm' for mouse, and 'a' for all markers.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return list of Seurat object that have the assigned clusters, and top 3 proposed cell types.
#'
#' @export
#'
sccca <- function(sobj, assay = "RNA", cluster, marker, tissue, tt = "a", cond, m_t = 0.9, c_t = 0.7, test = "p", org = 'a'){

  cell_clusters_names <- select(sobj@meta.data, cluster)
  #Load the database
  markers_db <- read.csv(marker)
  if (org == 'h'){
    markers_db <- markers_db[markers_db$species == "Human",]
  }else if (org == 'm')
  {
    markers_db <- markers_db[markers_db$species == "Mouse",]
  }
  if (tt == "n"){
    markers_db <- markers_db[!grepl(x = markers_db$celltype, pattern = "Leukemia|Cancer"),]
  }else if (tt == "c"){
    markers_db <- markers_db[grepl(x = markers_db$celltype, pattern = "Leukemia|Cancer"),]
  }
  # Get the organism tissue
  markers_db_tissue <- markers_db[grepl(tissue, markers_db$organ),]
  # Process the markers
  markers_list <- process_markers(markers_db_tissue)
  # Iterate over the cell types and check the markers
  cell_clusters <- unlist(unique(sobj@meta.data[cluster]))

  if(!is.null(cond)){
    cor_mat <- calculate_cor_mat(sobj, condition = cond, clusters = cluster)
  }else{
    cor_mat <- calculate_cor_mat(sobj, condition = NULL, clusters = cluster)
  }
  intersection_size <- NULL
  proposed_cell_types <- data.frame(matrix(nrow = 1, ncol = 4))
  colnames(proposed_cell_types) <- c("p.value", "intersection_size", "cluster", "cell_type")

  proposed_cell_types <- do.call(rbind, lapply(cell_clusters, process_clus,
                                               sobj = sobj, clus = cell_clusters_names, assay = assay,
                                               markers = markers_list, cor_m = cor_mat, m_t = m_t, c_t = c_t, test = test))
  proposed_cell_types_filt <- proposed_cell_types %>%
                                   group_by(cluster) %>%
                                   top_n(wt = intersection_size, n = 1)
  proposed_cell_types_filt <- proposed_cell_types_filt[!duplicated(proposed_cell_types_filt$cluster),]
  sobj@meta.data$cell_id <- rownames(sobj@meta.data)
  meta.data <- merge(y=sobj@meta.data, x=proposed_cell_types_filt,
                    by.y = cluster, by.x = "cluster")
  meta.data <- meta.data[order(meta.data$cell_id),]
  sobj@meta.data <- sobj@meta.data[order(rownames(sobj@meta.data)),]
  sobj@meta.data$cluster <- meta.data$cluster
  sobj@meta.data$cell_type <- meta.data$cell_type
  sobj@meta.data$p.value <- meta.data$p.value
  sobj@meta.data$intersection_size <- meta.data$intersection_size
  return(list(sobj, proposed_cell_types))
}
