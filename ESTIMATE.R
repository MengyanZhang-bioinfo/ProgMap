# ============================================================
# ESTIMATE Scoring Analysis
# ============================================================

cat("Performing ESTIMATE scoring analysis...\n")

# Prepare expression data for ESTIMATE
prepare_estimate_input <- function(expression_data, output_dir) {
  cat("Preparing ESTIMATE input data...\n")
  
  # Convert to log2 CPM
  dat <- log2(edgeR::cpm(expression_data) + 1)
  
  # Write input file
  write.table(dat, 
              file = file.path(output_dir, "CAN_estimate_input.txt"), 
              sep = '\t', quote = FALSE)
  
  return(dat)
}

estimate_input <- prepare_estimate_input(combined_expression, OUTPUT_DIR)

# 12.2 Run ESTIMATE analysis
run_estimate_analysis <- function(input_file, output_dir) {
  cat("Running ESTIMATE analysis...\n")
  
  # Create GCT format file
  in_gct_file <- file.path(output_dir, "ESTIMATE_CAN_input.gct")
  outputGCT(input_file, in_gct_file)
  
  # Filter common genes
  out_gene_file <- file.path(output_dir, "ESTIMATE_CAN_gene.gct")
  filterCommonGenes(input.f = input_file,
                    output.f = out_gene_file,
                    id = "GeneSymbol")
  
  # Calculate scores
  out_score_file <- file.path(output_dir, "ESTIMATE_score.gct")
  estimateScore(in_gct_file,
                out_score_file,
                platform = "illumina")
  
  # Read and process results
  estimate_scores <- read.table(out_score_file,
                                skip = 2,
                                header = TRUE,
                                row.names = 1)
  
  estimate_scores <- as.data.frame(t(estimate_scores[, 2:ncol(estimate_scores)]))
  
  # Save results
  write.table(estimate_scores,
              file = file.path(output_dir, "CAN_estimate_score.txt"),
              sep = '\t', quote = FALSE)
  
  return(estimate_scores)
}

estimate_scores <- run_estimate_analysis(file.path(OUTPUT_DIR, "CAN_estimate_input.txt"), 
                                         OUTPUT_DIR)

#Prepare ESTIMATE data for visualization
prepare_estimate_visualization_data <- function(estimate_scores, cluster_data) {
  cat("Preparing ESTIMATE data for visualization...\n")
  
  estimate_scores$sample <- substr(rownames(estimate_scores), 1, 12)
  
  # Match with cluster data
  sample_cluster_map <- data.frame(
    sample = substr(gsub('\\.', '-', rownames(cluster_df)), 1, 12),
    cluster = as.character(cluster_df$final_clusters)
  )
  
  estimate_merged <- merge(estimate_scores, sample_cluster_map, by = "sample")
  
  return(estimate_merged)
}

estimate_visualization_data <- prepare_estimate_visualization_data(estimate_scores, cluster_df)

# 12.4 Visualize ESTIMATE scores
visualize_estimate_scores <- function(estimate_data, output_file) {
  cat("Creating ESTIMATE score visualizations...\n")
  
  # Create individual plots
  stromal_plot <- ggplot(estimate_data, aes(x = cluster, y = StromalScore, fill = cluster)) +
    geom_boxplot(position = position_dodge(0.8)) +
    scale_fill_manual(values = CLUSTER_COLORS) +
    labs(x = "", y = 'Stromal Score') +
    stat_compare_means() +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75),
          legend.position = 'none')
  
  immune_plot <- ggplot(estimate_data, aes(x = cluster, y = ImmuneScore, fill = cluster)) +
    geom_boxplot(position = position_dodge(0.8)) +
    scale_fill_manual(values = CLUSTER_COLORS) +
    labs(x = "", y = 'Immune Score') +
    stat_compare_means() +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75),
          legend.position = 'none')
  
  estimate_plot <- ggplot(estimate_data, aes(x = cluster, y = ESTIMATEScore, fill = cluster)) +
    geom_boxplot(position = position_dodge(0.8)) +
    scale_fill_manual(values = CLUSTER_COLORS) +
    labs(x = "", y = 'ESTIMATE Score') +
    stat_compare_means() +
    theme_bw(base_size = 16)
  
  # Combine plots
  combined_plot <- stromal_plot + immune_plot + estimate_plot + 
    plot_layout(ncol = 3)
  
  # Save combined plot
  ggsave(output_file, combined_plot, width = 12, height = 5)
  
  return(list(
    stromal = stromal_plot,
    immune = immune_plot,
    estimate = estimate_plot,
    combined = combined_plot
  ))
}

estimate_plots <- visualize_estimate_scores(estimate_visualization_data,
                                            file.path(OUTPUT_DIR, "CAN_estimate_scores_plot.pdf"))

cat("ESTIMATE scoring analysis completed.\n")