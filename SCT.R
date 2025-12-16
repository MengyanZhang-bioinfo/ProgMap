##################################################
##################################################

analyze_spatial_transcriptomics <- function(data_dir, sample_id = "GSM") {
  # Create sample directory
  sample_dir <- file.path(data_dir, sample_id)
  dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
  setwd(sample_dir)
  
  # Load spatial data
  spatial_image <- Read10X_Image(
    file.path(sample_dir, "Spatial"),
    image.name = "tissue_hires_image.png"
  )
  
  spatial_obj <- Load10X_Spatial(
    data.dir = sample_dir,
    filename = paste0(sample_id, "_s1_filtered_feature_bc_matrix.h5"),
    assay = "Spatial",
    slice = "hnsc",
    image = spatial_image
  )
  
  # Quality control
  mt_genes <- grep(pattern = "^MT-", x = rownames(spatial_obj), value = TRUE)
  spatial_obj$percent.mito <- (
    Matrix::colSums(spatial_obj@assays[["Spatial"]]$counts[mt_genes, ]) /
      Matrix::colSums(spatial_obj@assays[["Spatial"]]$counts)
  ) * 100
  
  # Filter genes and spots
  genes_to_keep <- setdiff(
    names(which(Matrix::rowSums(spatial_obj@assays[["Spatial"]]$counts) > 5)),
    mt_genes
  )
  
  spatial_obj <- subset(
    spatial_obj,
    features = genes_to_keep,
    subset = nFeature_Spatial > 300 & percent.mito < 30
  )
  
  # Normalization
  spatial_obj <- SCTransform(spatial_obj, assay = "Spatial", verbose = FALSE)
  
  # Dimensionality reduction and clustering
  spatial_obj <- RunPCA(spatial_obj, assay = "SCT", verbose = FALSE)
  spatial_obj <- FindNeighbors(spatial_obj, reduction = "pca", dims = 1:30)
  spatial_obj <- FindClusters(spatial_obj, verbose = FALSE)
  spatial_obj <- RunUMAP(spatial_obj, reduction = "pca", dims = 1:30)
  
  return(spatial_obj)
}

visualize_spatial_data <- function(spatial_obj, marker_genes, output_dir) {
  # UMAP and spatial plots
  p_umap <- DimPlot(spatial_obj, reduction = "umap", label = TRUE)
  p_spatial <- SpatialDimPlot(spatial_obj, label = TRUE, label.size = 3, pt.size.factor = 4)
  
  # Marker gene visualization
  p_dotplot <- DotPlot(spatial_obj, features = marker_genes) + coord_flip()
  p_violin <- VlnPlot(spatial_obj, features = marker_genes, pt.size = 0, ncol = 5)
  
  # Save plots
  ggsave(file.path(output_dir, "spatial_umap.pdf"), p_umap + p_spatial, 
         width = 10, height = 10)
  ggsave(file.path(output_dir, "marker_dotplot.pdf"), p_dotplot, 
         width = 10, height = 10)
  ggsave(file.path(output_dir, "marker_violin.pdf"), p_violin, 
         width = 10, height = 10)
  
  return(list(
    umap_plot = p_umap,
    spatial_plot = p_spatial,
    dotplot = p_dotplot,
    violin_plot = p_violin
  ))
}

##################################################
##################################################

main_analysis_pipeline <- function() {
  # Load and preprocess data
  cat("Loading and preprocessing data...\n")
  seurat_obj <- load_scRNA_data()
  seurat_obj <- normalize_and_scale_data(seurat_obj)
  
  # Clustering
  cat("Performing clustering...\n")
  seurat_obj <- perform_clustering(seurat_obj)
  cluster_plots <- visualize_clusters(seurat_obj)
  
  # Cell type annotation
  cat("Annotating cell types...\n")
  annotated_obj <- manual_celltype_annotation(
    seurat_obj, 
    "G:/pan_result/gpu_result_tcga/CAN/marker_celltype.csv"
  )
  
  # InferCNV analysis
  cat("Running InferCNV...\n")
  infercnv_results <- run_infercnv(annotated_obj$obj)
  malignant_cells <- identify_malignant_cells("inferCNV_results")
  
  #Pseudotime analysis
  cat("Performing pseudotime analysis...\n")
  expressed_genes <- read.csv("./featuresweightme.csv", 
                             row.names = 1)
  expressed_genes <- colnames(expressed_genes)
  
  cds <- run_pseudotime_analysis(annotated_obj$obj, expressed_genes)
  
  # Define cell type colors
  cell_type_colors <- c(
    "B_cell" = "#59A14F",
    "Endothelium" = "#4E79A7",
    "Epithelium" = "#F28E2B",
    "Fibroblasts" = "#E15759",
    "Macrophage" = '#8DD3C7',
    "T_cell" = "#815c94",
    "Mast" = "#C8BF2C",
    "NK" = "#8c310a",
    "Plasmocyte" = '#BEBADA'
  )
  
  pseudotime_plots <- visualize_pseudotime(cds, color_palette, cell_type_colors)
  
  # Differential expression analysis
  cat("Performing differential expression analysis...\n")
  gene_data <- read.csv("./featuresweightme.csv", 
                       row.names = 1)
  cluster_labels <- read.table("./CAN_cluster.txt")
  
  deg_results <- perform_differential_expression(gene_data, cluster_labels)
  
  # AUCell analysis
  cat(" Running AUCell analysis...\n")
  gene_sets <- list(
    "Full_gene_set" = expressed_genes
  )
  
  aucell_results <- run_aucell_analysis(cds, gene_sets)
  
  
  # Spatial transcriptomics analysis
  cat("Analyzing spatial transcriptomics data...\n")
  spatial_obj <- analyze_spatial_transcriptomics(
    "./SCT"
  )
  
  marker_genes <- read.csv("./single_cell_marker_genes.csv")$marker
  spatial_plots <- visualize_spatial_data(
    spatial_obj, 
    marker_genes,
    "./GSM"
  )
  
  cat("Analysis pipeline completed successfully!\n")
  
  return(list(
    seurat_obj = seurat_obj,
    annotated_obj = annotated_obj,
    infercnv_results = infercnv_results,
    cds = cds,
    deg_results = deg_results,
    aucell_results = aucell_results,
    correlation_results = correlation_results,
    spatial_obj = spatial_obj
  ))
}

##################################################
##################################################

# Run the complete analysis pipeline
results <- main_analysis_pipeline()

# Optional: Save all results
save(results, file = "CAN_complete_analysis_results.RData")