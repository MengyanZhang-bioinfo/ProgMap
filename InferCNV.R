##################################################
##################################################

run_infercnv <- function(seurat_obj, output_dir = "inferCNV_results") {
  # Set up directories
  dir.create(output_dir, showWarnings = FALSE)
  setwd(output_dir)
  
  # Extract epithelial cells
  seurat_obj$cell_type <- seurat_obj@active.ident
  epithelial_cells <- subset(seurat_obj, idents = c('Epithelium'))
  
  # Prepare input files for inferCNV
  # 1. Expression matrix
  expr_matrix <- as.data.frame(
    GetAssayData(subset(seurat_obj, cells = colnames(epithelial_cells)))
  )
  
  # 2. Cell annotation file
  cell_annotations <- data.frame(
    cell_id = colnames(expr_matrix),
    cell_type = epithelial_cells@meta.data$cell_type
  )
  
  # 3. Gene information file
  gene_info <- annoGene(rownames(expr_matrix), "SYMBOL", 'human')
  gene_info <- gene_info[with(gene_info, order(chr, start)), c(1, 4:6)]
  gene_info <- gene_info[!duplicated(gene_info[, 1]), ]
  
  # Filter and sort expression matrix
  filtered_expr <- expr_matrix[rownames(expr_matrix) %in% gene_info[, 1], ]
  filtered_expr <- filtered_expr[match(gene_info[, 1], rownames(filtered_expr)), ]
  
  # Save input files
  write.table(filtered_expr, file = "expression_matrix.txt", sep = '\t', quote = FALSE)
  write.table(cell_annotations, file = "cell_annotations.txt", 
              sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(gene_info, file = "gene_info.txt", 
              sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # Run inferCNV
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = as.matrix(filtered_expr),
    annotations_file = "cell_annotations.txt",
    delim = "\t",
    gene_order_file = "gene_info.txt",
    ref_group_names = NULL
  )
  
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = output_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    num_threads = 8
  )
  
  return(infercnv_obj)
}

identify_malignant_cells <- function(infercnv_dir) {
  # Read CNV predictions
  cnv_regions <- read.table(
    file.path(infercnv_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat"),
    header = TRUE
  )
  
  cell_groupings <- read.table(
    file.path(infercnv_dir, "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"),
    header = TRUE
  )
  
  # Identify malignant cells
  malignant_cells <- cell_groupings[
    cell_groupings$cell_group_name %in% cnv_regions$cell_group,
  ]
  
  # Save malignant cell list
  malignant_cell_ids <- malignant_cells$cell
  write.table(malignant_cell_ids, "CAN_EPI_cells_E.txt")
  
  return(malignant_cell_ids)
}