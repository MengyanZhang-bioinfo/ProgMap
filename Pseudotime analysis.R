##################################################
##################################################

run_pseudotime_analysis <- function(seurat_obj, expressed_genes, output_dir = "pseudotime_analysis") {
  dir.create(output_dir, showWarnings = FALSE)
  setwd(output_dir)
  
  # Create monocle object
  expr_matrix <- as(as.matrix(seurat_obj@assays[["RNA"]]@layers[["counts"]]), 'sparseMatrix')
  p_data <- seurat_obj@meta.data
  f_data <- data.frame(
    gene_short_name = rownames(seurat_obj),
    row.names = rownames(seurat_obj)
  )
  
  pd <- new('AnnotatedDataFrame', data = p_data)
  fd <- new('AnnotatedDataFrame', data = f_data)
  
  cds <- newCellDataSet(
    expr_matrix,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0.1,
    expressionFamily = negbinomial.size()
  )
  
  # Process monocle object
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- detectGenes(cds, min_expr = 0.1)
  
  # Set ordering genes
  cds <- setOrderingFilter(cds, expressed_genes)
  
  # Reduce dimensions and order cells
  cds <- reduceDimension(cds, max_components = 5, method = 'DDRTree')
  cds <- orderCells(cds)
  
  # Save cds object
  save(cds, file = "cds.RData")
  
  return(cds)
}

visualize_pseudotime <- function(cds, color_palette, cell_type_colors) {
  # Create pseudotime plots
  p_pseudotime <- plot_cell_trajectory(cds, color_by = "Pseudotime", 
                                       size = 1, show_backbone = TRUE)
  
  p_celltype <- plot_cell_trajectory(cds, color_by = "cell_type", 
                                     size = 1, show_backbone = TRUE) +
    scale_color_manual(values = cell_type_colors) +
    theme(legend.position = "right")
  
  # Create density plot of pseudotime by cell type
  df <- pData(cds)
  p_density <- ggplot(df, aes(Pseudotime, colour = cell_type, fill = cell_type)) +
    geom_density(bw = 0.5, size = 1, alpha = 0.5) +
    theme_classic2()
  
  # Save plots
  pdf("pseudotime_analysis.pdf", width = 12, height = 10)
  print(p_pseudotime)
  print(p_celltype)
  print(p_density)
  dev.off()
  
  return(list(
    pseudotime_plot = p_pseudotime,
    celltype_plot = p_celltype,
    density_plot = p_density
  ))
}

##################################################

##################################################

run_aucell_analysis <- function(cds, gene_sets) {
  # Extract expression matrix
  expr_matrix <- exprs(cds)
  
  # Convert to matrix if needed
  if (class(expr_matrix) == "dgCMatrix") {
    expr_matrix <- as.matrix(expr_matrix)
  }
  
  # Build rankings
  set.seed(123)
  cells_rankings <- AUCell_buildRankings(
    expr_matrix,
    nCores = 1,
    plotStats = TRUE
  )
  
  # Calculate AUC values
  cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)
  
  # Get AUC matrix
  auc_matrix <- getAUC(cells_AUC)
  
  # Normalize AUC values to 0-1 range
  normalize_0_1 <- function(x) {
    if (max(x) == min(x)) return(rep(0.5, length(x)))
    (x - min(x)) / (max(x) - min(x))
  }
  
  auc_matrix_normalized <- t(apply(auc_matrix, 1, normalize_0_1))
  rownames(auc_matrix_normalized) <- rownames(auc_matrix)
  colnames(auc_matrix_normalized) <- colnames(auc_matrix)
  
  # Add AUC scores to cds object
  cell_names <- colnames(cds)
  
  for (gene_set_name in rownames(auc_matrix_normalized)) {
    # Normalized scores
    auc_scores_norm <- auc_matrix_normalized[gene_set_name, cell_names]
    pData(cds)[[paste0("AUC_", gene_set_name, "_norm")]] <- auc_scores_norm
    
    # Original scores
    auc_scores_orig <- auc_matrix[gene_set_name, cell_names]
    pData(cds)[[paste0("AUC_", gene_set_name, "_orig")]] <- auc_scores_orig
  }
  
  return(list(
    cds = cds,
    auc_matrix = auc_matrix,
    auc_matrix_normalized = auc_matrix_normalized
  ))
}



