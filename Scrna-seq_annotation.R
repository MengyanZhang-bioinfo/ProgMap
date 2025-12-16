##################################################
##################################################
rm(list = ls())
# Install required packages if not installed
# install.packages(c("Seurat", "cowplot", "dplyr", "data.table", "ggthemes", "ggplot2", 
#                    "clustree", "patchwork", "tidyr", "ggpubr", "viridis", "gridExtra"))
# BiocManager::install(c("SingleR", "celldex", "scater", "monocle", "infercnv", 
#                       "AUCell", "GSEABase", "BayesPrism", "AnnoProbe", "rjags"))

# Load libraries
library(Seurat)
library(cowplot)
library(dplyr)
library(data.table)
library(ggthemes)
library(ggplot2)
library(clustree)
library(patchwork)
library(tidyr)
library(ggpubr)
library(viridis)
library(gridExtra)

# Bioconductor libraries
library(SingleR)
library(celldex)
library(scater)
library(monocle)
library(infercnv)
library(AUCell)
library(GSEABase)
library(BayesPrism)
library(AnnoProbe)

# Set working directory
setwd("G:/pan_result/gpu_result_tcga/CAN/scRNA_new")

##################################################
##################################################

load_scRNA_data <- function() {
  # Read count matrix
  meta <- read.table(gzfile("_UMI_counts.txt.gz"), header = TRUE, sep = "\t")
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = meta,
    min.features = 200,
    min.cells = 3,
    project = "GSE"
  )
  
  # Quality control
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > 200 & nCount_RNA > 2500 & percent.mt < 20)
  
  return(seurat_obj)
}

normalize_and_scale_data <- function(seurat_obj) {
  # Normalization
  seurat_obj <- NormalizeData(seurat_obj, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000)
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, 
                                     selection.method = "vst", 
                                     nfeatures = 2000)
  
  # Scale data
  seurat_obj <- ScaleData(seurat_obj)
  
  # PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), dims = 1:30)
  
  return(seurat_obj)
}

##################################################
##################################################

perform_clustering <- function(seurat_obj, pc_num = 1:50, resolution = 0.6) {
  # Run dimensionality reduction
  seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = pc.num) %>%
    RunUMAP(reduction = "pca", dims = pc.num) %>%
    FindNeighbors(reduction = "pca", dims = pc.num)
  
  # Find clusters
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  return(seurat_obj)
}

visualize_clusters <- function(seurat_obj, output_dir = "harmony_0.4_20") {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE)
  setwd(output_dir)
  
  # Define color palette
  color_palette <- c("#59A14F", "#4E79A7", "#F28E2B", "#E15759", "#4E80AB",
                     "#FFCC33", "#C8BF2C", "#8c310a", "#228B22", "#815c94",
                     "#3692a8", "#c95043", "#b5574d", "#548faf", "#e97e33",
                     "#76a695", '#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', 
                     '#80B1D3', '#FDB462')
  
  # Create UMAP plots
  p_umap_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, raster = FALSE)
  p_umap_samples <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident",
                           split.by = "orig.ident", ncol = 3, cols = color_palette, 
                           raster = FALSE) + NoLegend()
  
  # Create tSNE plots
  p_tsne_clusters <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE, raster = FALSE)
  p_tsne_samples <- DimPlot(seurat_obj, reduction = "tsne", group.by = "orig.ident",
                           split.by = "orig.ident", ncol = 3, cols = color_palette,
                           raster = FALSE) + NoLegend()
  
  # Save plots
  ggsave("UMAP_Samples_harmony.pdf", p_umap_samples, width = 9, height = 9)
  ggsave("UMAP_cluster_harmony.pdf", p_umap_clusters, width = 9, height = 9)
  ggsave("TSNE_Samples_harmony.pdf", p_tsne_samples, width = 9, height = 9)
  ggsave("TSNE_cluster_harmony.pdf", p_tsne_clusters, width = 9, height = 9)
  
  # Save Seurat object
  saveRDS(seurat_obj, "Seurat.rds")
  
  # Save cluster information
  dir.create("Markers", showWarnings = FALSE)
  cluster_info <- data.frame(Cluster = seurat_obj@active.ident)
  write.csv(cluster_info, 'Markers/Cluster.csv', quote = FALSE)
  
  return(list(
    umap_clusters = p_umap_clusters,
    umap_samples = p_umap_samples,
    tsne_clusters = p_tsne_clusters,
    tsne_samples = p_tsne_samples
  ))
}

##################################################
##################################################

manual_celltype_annotation <- function(seurat_obj, marker_file) {
  # Load manual annotation file
  manual_celltype <- read.csv(marker_file)
  
  # Rename clusters
  new.cluster.ids <- manual_celltype$celltype
  names(new.cluster.ids) <- levels(seurat_obj)
  seurat_annotated <- RenameIdents(seurat_obj, new.cluster.ids)
  
  # Create visualization plots
  p_umap_annotated <- DimPlot(seurat_annotated, reduction = "umap", 
                              label = TRUE, pt.size = 0.5) + NoLegend()
  p_tsne_annotated <- DimPlot(seurat_annotated, reduction = "tsne", 
                              label = TRUE, pt.size = 0.5) + NoLegend()
  
  # Save annotated object
  save(seurat_annotated, file = "GSE.seu.obj.Rdata")
  
  return(list(
    obj = seurat_annotated,
    umap_plot = p_umap_annotated,
    tsne_plot = p_tsne_annotated
  ))
}

singleR_annotation <- function(seurat_obj, ref_database_path) {
  # Load reference database
  ref_data <- get(load(ref_database_path))
  
  # Prepare data for SingleR
  seurat_obj$cell.type <- Idents(seurat_obj)
  sce_obj <- as.SingleCellExperiment(seurat_obj)
  
  # Run SingleR
  singler_result <- SingleR(
    test = sce_obj,
    ref = ref_data,
    labels = ref_data$label.main,
    method = "cluster",
    cluster = sce_obj$cell.type
  )
  
  # Extract annotation
  annotation_df <- singler_result %>% 
    dplyr::tbl_df() %>% 
    dplyr::select(cluster, labels)
  
  # Apply annotation
  new.cluster.ids <- annotation_df$labels
  names(new.cluster.ids) <- annotation_df$cluster
  seurat_annotated <- RenameIdents(seurat_obj, new.cluster.ids)
  
  # Save results
  dir.create("singleR", showWarnings = FALSE)
  write.table(annotation_df, file = "singleR/cluster_after.xls", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(list(
    obj = seurat_annotated,
    annotation = annotation_df,
    singler_result = singler_result
  ))
}
