# ============================================================
# TIDE Analysis
# ============================================================

cat("Performing TIDE analysis...\n")

#Load TIDE results
tide_results <- read.csv(TIDE_FILE, header = TRUE)

cluster_df <- as.data.frame(final_clusters)
cluster_tide_data <- data.frame(
  Patient = rownames(cluster_df),
  Group = cluster_df$final_clusters
)

# Merge with TIDE results
tide_merged <- left_join(tide_results, cluster_tide_data, by = "Patient")

#Statistical test
tide_test <- wilcox.test(as.numeric(tide_merged$TIDE) ~ tide_merged$Group)
tide_pvalue <- tide_test$p.value

cat(sprintf("TIDE analysis p-value: %.4f\n", tide_pvalue))

# Create TIDE visualization
tide_plot_data <- data.frame(
  TIDE = as.numeric(tide_merged$TIDE),
  Cluster = as.character(tide_merged$Group)
)

# Create advanced violin-box plot
create_tide_plot <- function(plot_data, output_file) {
  # Base plot
  base_plot <- ggplot(plot_data, aes(x = Cluster, y = TIDE, fill = Cluster, color = Cluster)) +
    scale_fill_manual(values = CLUSTER_COLORS) +
    scale_colour_manual(values = CLUSTER_COLORS)
  
  # Add half violin
  plot_violin <- base_plot + geom_half_violin(
    position = position_nudge(x = 0.1, y = 0),
    side = 'R',
    adjust = 1.2,
    trim = FALSE,
    color = NA,
    alpha = 0.6
  )
  
  # Add half points
  plot_points <- plot_violin + geom_half_point(
    position = position_nudge(x = -0.35, y = 0),
    size = 3,
    shape = 19,
    range_scale = 0.5,
    alpha = 0.6
  )
  
  # Add boxplot
  plot_box <- plot_points + geom_boxplot(
    outlier.shape = NA,
    width = 0.1,
    alpha = 0.8
  )
  
  # Add statistical comparison
  plot_stats <- plot_box + stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.x = 1.5,
    size = 5
  )
  
  # Final theme customization
  final_plot <- plot_stats + theme_light() + 
    theme(panel.grid = element_blank()) +
    theme(
      axis.ticks.length = unit(-0.25, "cm"),
      axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
      axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
      panel.grid.major = element_line(linetype = "dotted", color = "lightgray", size = 0.5),
      text = element_text(size = 20)
    ) +
    coord_cartesian(ylim = c(-4, 4))
  
  # Save plot
  pdf(output_file, width = 7, height = 7)
  print(final_plot)
  dev.off()
  
}

# Create and save TIDE plot
cluster_tide_plot <- create_tide_plot(tide_plot_data, 
                                      file.path(OUTPUT_DIR, "CAN_cluster_tide_plot.pdf"))

cat("TIDE analysis completed.\n")