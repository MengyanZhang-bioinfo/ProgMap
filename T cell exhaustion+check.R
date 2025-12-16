# ============================================================
# T-cell Exhaustion Analysis
# ============================================================

cat("Performing T-cell exhaustion analysis...\n")

# Load and prepare checkpoint genes
load_checkpoint_genes <- function(checkpoint_file) {
  cat("Loading checkpoint genes...\n")
  
  checkpoint_data <- read.csv(checkpoint_file)
  return(checkpoint_data)
}

checkpoint_genes <- load_checkpoint_genes("./checkpoint_genes.csv")

#Analyze checkpoint gene expression
analyze_checkpoint_expression <- function(expression_data, checkpoint_data, cluster_samples) {
  cat("Analyzing checkpoint gene expression...\n")
  
  # Get checkpoint genes present in expression data
  checkpoint_in_expr <- checkpoint_data$Checkpoint.gene[checkpoint_data$Checkpoint.gene %in% rownames(expression_data)]
  checkpoint_expr <- expression_data[checkpoint_in_expr, ]
  checkpoint_expr <- na.omit(checkpoint_expr)
  
  # Separate by cluster
  cluster1_check <- checkpoint_expr[, colnames(checkpoint_expr) %in% cluster_samples[[1]]]
  cluster2_check <- checkpoint_expr[, colnames(checkpoint_expr) %in% cluster_samples[[2]]]
  
  # Combine for analysis
  combined_check <- cbind(cluster1_check, cluster2_check)
  
  # Statistical analysis
  p_values <- numeric(nrow(combined_check))
  cluster_groups <- rep(c(1, 2), c(ncol(cluster1_check), ncol(cluster2_check)))
  
  for (j in 1:nrow(combined_check)) {
    test_result <- kruskal.test(as.numeric(combined_check[j, ]) ~ cluster_groups)
    p_values[j] <- test_result$p.value
  }
  
  # Identify significant checkpoint genes
  significant_indices <- which(p_values < 0.05)
  significant_genes <- rownames(combined_check)[significant_indices]
  
  cat(sprintf("Found %d significant checkpoint genes (p < 0.05)\n", length(significant_genes)))
  
  # Prepare detailed results
  significant_details <- checkpoint_data[checkpoint_data$Checkpoint.gene %in% significant_genes, ]
  significant_details <- significant_details[order(significant_details$Activator.A..Inhibitor.I.), ]
  
  # Calculate mean expression
  cluster1_mean <- apply(checkpoint_expr[significant_details$Checkpoint.gene, colnames(cluster1_check)], 1, mean)
  cluster2_mean <- apply(checkpoint_expr[significant_details$Checkpoint.gene, colnames(cluster2_check)], 1, mean)
  
  return(list(
    all_checkpoints = checkpoint_expr,
    significant_genes = significant_genes,
    significant_details = significant_details,
    p_values = p_values,
    cluster1_mean = cluster1_mean,
    cluster2_mean = cluster2_mean,
    cluster1_expression = cluster1_check,
    cluster2_expression = cluster2_check
  ))
}

checkpoint_analysis <- analyze_checkpoint_expression(combined_expression, checkpoint_genes, cluster_samples_list)

# 14.3 Save checkpoint analysis results
save_checkpoint_results <- function(checkpoint_analysis, output_dir) {
  cat("Saving checkpoint analysis results...\n")
  
  # Save significant genes with p-values
  sig_genes_p <- cbind(checkpoint_analysis$significant_genes,
                       checkpoint_analysis$p_values[which(checkpoint_analysis$p_values < 0.05)])
  
  write.table(sig_genes_p,
              file.path(output_dir, "CAN_checkpoint_significant_genes.txt"),
              row.names = FALSE, col.names = c("Gene", "p_value"))
  
  # Save detailed information
  write.table(checkpoint_analysis$significant_details,
              file.path(output_dir, "CAN_checkpoint_detailed_info.txt"),
              row.names = FALSE)
}

save_checkpoint_results(checkpoint_analysis, OUTPUT_DIR)

# Create radar plot for checkpoint genes
create_checkpoint_radar_plot <- function(cluster1_means, cluster2_means, gene_names, output_file) {
  cat("Creating checkpoint gene radar plot...\n")
  
  # Prepare data for radar plot
  radar_data <- rbind(cluster1_means, cluster2_means)
  rownames(radar_data) <- c("cluster1", "cluster2")
  
  # Convert to dataframe with cluster labels
  radar_df <- data.frame(
    cluster = c("cluster1", "cluster2"),
    as.data.frame(t(radar_data))
  )
  
  # Convert to numeric
  radar_df[, 2:ncol(radar_df)] <- apply(radar_df[, 2:ncol(radar_df)], 2, as.numeric)
  
  # Determine scale limits
  min_value <- min(radar_df[, 2:ncol(radar_df)], na.rm = TRUE)
  max_value <- max(radar_df[, 2:ncol(radar_df)], na.rm = TRUE)
  mid_value <- (min_value + max_value) / 2
  
  # Create radar plot
  radar_plot <- ggradar(plot.data = radar_df,
                        grid.min = floor(min_value * 10) / 10,
                        grid.mid = round(mid_value, 2),
                        grid.max = ceiling(max_value * 10) / 10,
                        values.radar = c(round(min_value, 2), round(mid_value, 2), round(max_value, 2)),
                        group.colours = CLUSTER_COLORS,
                        axis.label.size = 4)
  
  ggsave(output_file, radar_plot, width = 10, height = 10)
  
  return(radar_plot)
}

# Create radar plot if significant genes found
if (length(checkpoint_analysis$significant_genes) > 0) {
  radar_plot <- create_checkpoint_radar_plot(
    checkpoint_analysis$cluster1_mean[checkpoint_analysis$significant_details$Checkpoint.gene],
    checkpoint_analysis$cluster2_mean[checkpoint_analysis$significant_details$Checkpoint.gene],
    checkpoint_analysis$significant_details$Checkpoint.gene,
    file.path(OUTPUT_DIR, "CAN_checkpoint_radar_plot.pdf")
  )
}

# T-cell exhaustion scoring (TCellSI)
perform_tcell_exhaustion_analysis <- function(expression_data, cluster_samples) {
  cat("Performing T-cell exhaustion analysis...\n")
  
  # Check if TCellSI is available
  if (!require("TCellSI", quietly = TRUE)) {
    cat("Installing TCellSI package...\n")
    devtools::install_github("GuoBioinfoLab/TCellSI")
    library(TCellSI)
  }
  
  # Calculate T-cell exhaustion scores
  exhaustion_scores <- TCellSI::TCSS_Calculate(expression_data)
  
  # Separate scores by cluster
  cluster1_scores <- exhaustion_scores[, colnames(exhaustion_scores) %in% cluster_samples[[1]]]
  cluster2_scores <- exhaustion_scores[, colnames(exhaustion_scores) %in% cluster_samples[[2]]]
  
  # Focus on progenitor and terminal exhaustion
  key_scores <- rbind(cluster1_scores, cluster2_scores)
  key_scores <- key_scores[c("Progenitor_exhaustion", "Terminal_exhaustion"), ]
  
  # Statistical analysis
  p_values <- numeric(nrow(key_scores))
  score_groups <- rep(c(2, 1), c(ncol(cluster2_scores), ncol(cluster1_scores)))
  
  for (j in 1:nrow(key_scores)) {
    test_result <- kruskal.test(as.numeric(key_scores[j, ]) ~ score_groups)
    p_values[j] <- test_result$p.value
  }
  
  cat(sprintf("Progenitor exhaustion p-value: %.4f\n", p_values[1]))
  cat(sprintf("Terminal exhaustion p-value: %.4f\n", p_values[2]))
  
  return(list(
    all_scores = exhaustion_scores,
    key_scores = key_scores,
    cluster1_scores = cluster1_scores,
    cluster2_scores = cluster2_scores,
    p_values = p_values
  ))
}

exhaustion_analysis <- perform_tcell_exhaustion_analysis(combined_expression, cluster_samples_list)

# 14.6 Visualize T-cell exhaustion results
visualize_exhaustion_scores <- function(exhaustion_data, output_file) {
  cat("Creating T-cell exhaustion visualization...\n")
  
  # Prepare data for plotting
  plot_data <- data.frame()
  score_types <- rep(c("cluster2", "cluster1"), 
                     c(ncol(exhaustion_data$cluster2_scores), 
                       ncol(exhaustion_data$cluster1_scores)))
  
  for (i in 1:nrow(exhaustion_data$key_scores)) {
    scores <- as.numeric(exhaustion_data$key_scores[i, ])
    score_type <- rep(rownames(exhaustion_data$key_scores)[i], length(scores))
    cluster <- score_types
    plot_data <- rbind(plot_data, data.frame(scores, score_type, cluster))
  }
  
  # Remove outliers
  plot_data <- plot_data %>%
    filter(scores < quantile(scores, 0.95))
  
  # Create violin-box plot
  exhaustion_plot <- ggplot(plot_data, aes(x = factor(score_type), y = scores, fill = cluster)) +
    geom_violin(position = position_dodge(0.8), trim = TRUE, alpha = 0.6) +
    geom_boxplot(width = 0.2, position = position_dodge(0.8), color = "black", 
                 outlier.shape = NA, alpha = 0.5) +
    scale_fill_manual(values = CLUSTER_COLORS) +
    labs(x = "", y = 'Exhaustion Score') +
    stat_compare_means(method = "kruskal.test", label = "p.signif") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
  
  ggsave(output_file, exhaustion_plot, width = 10, height = 10)
  
  return(exhaustion_plot)
}

exhaustion_plot <- visualize_exhaustion_scores(exhaustion_analysis,
                                               file.path(OUTPUT_DIR, "CAN_tcell_exhaustion_plot.pdf"))

cat("T-cell exhaustion analysis completed.\n")