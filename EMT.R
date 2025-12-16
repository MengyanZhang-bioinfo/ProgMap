
#EMT analysis
perform_emt_analysis <- function(expression_data, cluster_samples) {
  cat("Performing EMT analysis...\n")
  
  # Load EMT gene sets
  emt_file <- "./emt_markers.xlsx"
  emt_data <- read.xlsx(emt_file)
  
  # Create gene set list
  emt_gene_list <- split(as.matrix(emt_data)[, 1], emt_data[, 2])
  
  # Perform GSVA
  emt_gsva_results <- perform_immune_gsva(expression_data, cluster_samples, emt_gene_list)
  
  # Differential analysis
  emt_analysis <- analyze_differential_immune(emt_gsva_results)
  
  return(list(
    gsva_results = emt_gsva_results,
    analysis = emt_analysis
  ))
}

emt_results <- perform_emt_analysis(combined_expression, cluster_samples_list)

# Create EMT visualization
if (nrow(emt_results$analysis$significant) > 0) {
  emt_plot_data <- data.frame(
    scores = c(as.numeric(emt_results$analysis$significant[1, ])),
    type = "Cell Migration",
    cluster = rep(c("cluster2", "cluster1"), 
                  c(ncol(emt_results$gsva_results$cluster2), 
                    ncol(emt_results$gsva_results$cluster1)))
  )
  
  emt_plot <- ggviolin(emt_plot_data, 
                       x = "cluster", 
                       y = "scores",
                       fill = "cluster",
                       facet.by = "type",
                       alpha = 1,
                       width = 0.5,
                       ylab = "Normalized Expression",
                       xlab = FALSE,
                       add = "boxplot",
                       add.params = list(fill = "white", width = 0.1, linetype = 1)) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_fill_manual(values = CLUSTER_COLORS)
  
  ggsave(file.path(OUTPUT_DIR, "CAN_emt_analysis_plot.pdf"), 
         emt_plot, width = 10, height = 10)
}

cat("Immune analysis completed.\n")