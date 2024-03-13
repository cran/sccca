#' Process the cell markers names 
#'
#' This Function is used to return the cell markers names processed for the sctype approach
#'
#' @usage correct_gene_symbols(markers)
#'
#' @param markers list of unique cell markers.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return vector of genes names which overlap with the correlation matrix.
#'
#' @export
#'
correct_gene_symbols <- function(markers) {
  gene_symbols <- gsub(" ", "", markers)
  gene_symbols <- toupper(gene_symbols[gene_symbols != "NA" & gene_symbols != ""])
  gene_symbols <- sort(gene_symbols)
}
#' Process the database for the sctype approach 
#'
#' This Function is used to process the database that will be used for sctype approach
#'
#' @usage process_database(database_name = "sctype", org = 'a', tissue, tissue_type = 'n')
#'
#' @param database_name name of the database to be used that can be 'sctype' or 'UMD'.
#' 
#' @param org name of organism to be used that can be 'h' for human, 'm' for mouse, and 'a' for all markers. 
#' 
#' @param tissue specified tissue from which the data comes.
#' 
#' @param tissue_type tissue type whether 'a' for all types 'n' for normal tissues only or "c" for cancer tissues.
#'
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return vector of genes names which overlap with the correlation matrix.
#'
#' @export
process_database <- function(database_name = "sctype", org = 'a', tissue, tissue_type = 'n'){
  if (database_name == "sctype")
  {
    
    cell_markers = read.xlsx("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx")
    cell_markers = cell_markers[cell_markers$tissueType == tissue,] 
    cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
    cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    
    # correct gene symbols from the given DB (up-genes)
    cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
      
      markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
      markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
      markers_all = sort(markers_all)
      
      if(length(markers_all) > 0){
        suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
        paste0(markers_all, collapse=",")
      } else {
        ""
      }
    })
    
    # correct gene symbols from the given DB (down-genes)
    cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
      
      markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
      markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
      markers_all = sort(markers_all)
      
      if(length(markers_all) > 0){
        suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
        paste0(markers_all, collapse=",")
      } else {
        ""
      }
    })
    
    cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
    cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    
    gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
    gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
    
    return(list(gs_positive = gs, gs_negative = gs2))
    
  }else if (database_name == "UMD"){
    
    markers_db <- read.csv("https://gitlab.lcsb.uni.lu/mohamed.soudy/sccca/-/raw/main/Unified_markers_db.csv")
    if (org == 'h'){
      markers_db <- markers_db[markers_db$species == "Human",]
    }else if (org == 'm')
    {
      markers_db <- markers_db[markers_db$species == "Mouse",]
    }
    if (tissue_type == "n"){
      markers_db <- markers_db[!grepl(x = toupper(markers_db$celltype), pattern = "LEUKEMIA|CANCER|CARCINOMA"),]
    }else if (tissue_type == "c"){
      markers_db <- markers_db[grepl(x = toupper(markers_db$celltype), pattern = "LEUKEMIA|CANCER|CARCINOMA"),]
    }
    # Get the organism tissue
    markers_db_tissue <- markers_db[grepl(tissue, markers_db$organ),]
    # Process the markers
    gs_positive <- process_markers(markers_db_tissue)
    return(list(gs_positive=gs_positive))
  }#UMD 
  
}

#' Run the sctype approach as it's implemented by Ianevski, A., Giri, A.K. and Aittokallio, T.
#'
#' This Function is used to run the sctype approach with faster implementation
#'
#' @usage sctype(sobj,assay="RNA",tissue,tt="a",clus,org="a",scaled=T,database="sctype")
#'
#' @param sobj Seurat object.
#' 
#' @param assay assay to be used default is set to RNA.
#' 
#' @param tissue specified tissue from which the data comes.
#' 
#' @param tt tissue type whether 'a' for all types 'n' for normal tissues only or "c" for cancer tissues.
#' 
#' @param clus colname in the mata.data that have the cell cluster numbers.
#' 
#' @param org organism to be used that can be 'h' for human, 'm' for mouse, and 'a' for all markers.
#' 
#' @param scaled indicates whether the matrix is scaled (TRUE by default)
#'
#' @param database name of the database to be used that can be 'sctype' or 'UMD'
#' 
#' @author Mohmed Soudy \email{Mohamed.soudy@uni.lu} and Sohpie LE BARS \email{sophie.lebars@uni.lu} and Enrico Glaab \email{enrico.glaab@uni.lu}
#'
#' @return vector of genes names which overlap with the correlation matrix.
#'
#' @export
sctype <- function(sobj, assay = "RNA", tissue, tt = "a", clus, org = "a", scaled = T,database="sctype"){
  
  scRNAseqData <- as.matrix(sobj@assays[[assay]]$counts)
  
  #Get the database 
  markers <- process_database(database_name = database, org = org, tissue = tissue, tissue_type = tt)
  
  marker_positive <- markers$gs_positive
  if (database == "sctype"){
    marker_negative <- markers$gs_negative
  }else{
    marker_negative <- NULL
  }
  marker_stat = sort(table(unlist(marker_positive)), decreasing = T)
  marker_sensitivity = data.frame(score_marker_sensitivity = rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(marker_positive),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  #Marker genes that are found in our data 
  matched_markers <- marker_sensitivity[match_characters(marker_sensitivity$gene_, scRNAseqData),]
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(marker_positive)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  Z <- data.frame(Z)
  row_indices <- cell_markers_genes_score[,"gene_"]
  
  # Create a matrix of score_marker_sensitivity values to be multiplied
  multiplier <- cell_markers_genes_score[,"score_marker_sensitivity"]
  
  # Use matrix indexing to perform the multiplication
  Z[row_indices, ] <- Z[row_indices, ] * multiplier
  
  . = NULL
  # subselect only with marker genes
  Z = Z[unique(unlist(marker_positive)), ]
  # combine scores
  es <- names(marker_positive) %>%
    lapply(function(gss_) {
      sapply(1:ncol(Z), function(j) {
        gs_z <- Z[marker_positive[[gss_]], j]
        gz_2 <- Z[marker_negative[[gss_]], j] * -1
        sum_t1 <- sum(gs_z, na.rm = TRUE) / sqrt(sum(!is.na(gs_z)))
        sum_t2 <- sum(gz_2, na.rm = TRUE) / sqrt(sum(!is.na(gz_2)))
        if (is.na(sum_t2)) sum_t2 <- 0
        sum_t1 + sum_t2
      })
    }) %>%
    do.call(rbind, .)
  
  dimnames(es) = list(names(marker_positive), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  # Extract top cell types for each cluster
  cluster <- NULL
  scores <- NULL
  cL_resutls = do.call("rbind", lapply(unique(sobj@meta.data[[clus]]), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sobj@meta.data[sobj@meta.data[[clus]]==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sobj@meta.data[[clus]]==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  # Update sctype_scores$type where the condition is met
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
  
  sobj@meta.data$cell_type <- ""
  # Iterate over unique cluster values and update meta.data accordingly
  sapply(unique(sctype_scores$cluster), function(j) {
    cl_type <- sctype_scores[sctype_scores$cluster == j,]
    sobj@meta.data[sobj@meta.data[[clus]] == j, 'cell_type'] <<- as.character(cl_type$type[1])
  })
  return(sobj)
  
}
