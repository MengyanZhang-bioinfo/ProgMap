# ============================================================
# Mutation Analysis
# ============================================================

cat("Step 11: Performing mutation analysis...\n")

# 11.1 Set up mutation data directory
mutation_base_dir <- "./Simple_Nucleotide_Variation"
setwd(mutation_base_dir)

# 11.2 Load and process MAF files
process_maf_files <- function(mutation_dir) {
  cat("Loading MAF files...\n")
  
  # Find all MAF files
  maf_files <- dir(path = "Masked_Somatic_Mutation\\", 
                   pattern = "masked.maf.gz$", 
                   full.names = TRUE, 
                   recursive = TRUE)
  
  cat(sprintf("Found %d MAF files\n", length(maf_files)))
  
  # Read MAF files with error handling
  maf_data_list <- lapply(maf_files, function(x) {
    tryCatch({
      read.maf(x, isTCGA = TRUE)
    }, error = function(e) {
      message(paste("Error reading file:", x, "\n", e$message))
      return(NA)
    })
  })
  
  # Remove failed reads
  maf_data_clean <- maf_data_list[!sapply(maf_data_list, is.na)]
  
  # Merge MAF files
  if (length(maf_data_clean) > 0) {
    merged_maf <- merge_mafs(maf_data_clean)
    cat(sprintf("Successfully merged %d MAF files\n", length(maf_data_clean)))
    return(merged_maf)
  } else {
    stop("No valid MAF files found")
  }
}

# Process MAF files
merged_maf <- process_maf_files(mutation_base_dir)

# Save merged MAF
save(merged_maf, file = file.path(mutation_base_dir, "TCGA-CAN_SNP.Rdata"))

# 11.3 Prepare cluster data for mutation analysis
prepare_mutation_cluster_data <- function(cluster_df) {
  mutation_cluster_data <- data.frame(
    sample = substr(gsub('\\.', '-', rownames(cluster_df)), 1, 12),
    cluster = cluster_df$final_clusters
  )
  return(mutation_cluster_data)
}

mutation_cluster_data <- prepare_mutation_cluster_data(cluster_df)

# 11.4 Generate MAF summary plots
generate_maf_summaries <- function(maf_data, cluster_data, output_dir) {
  cat("Generating MAF summary plots...\n")
  
  # Overall MAF summary
  pdf(file.path(output_dir, "CAN_overall_maf_summary.pdf"))
  plotmafSummary(maf = maf_data, 
                 rmOutlier = TRUE, 
                 addStat = 'median', 
                 dashboard = TRUE)
  dev.off()
  
  # Subset by cluster
  cluster1_samples <- cluster_data$sample[cluster_data$cluster == "1"]
  cluster2_samples <- cluster_data$sample[cluster_data$cluster == "2"]
  
  maf_cluster1 <- subsetMaf(maf = maf_data, tsb = cluster1_samples)
  maf_cluster2 <- subsetMaf(maf = maf_data, tsb = cluster2_samples)
  
  # Cluster 1 summary
  pdf(file.path(output_dir, "CAN_cluster1_maf_summary.pdf"))
  plotmafSummary(maf = maf_cluster1, 
                 rmOutlier = TRUE, 
                 addStat = 'median', 
                 dashboard = TRUE)
  dev.off()
  
  # Cluster 2 summary
  pdf(file.path(output_dir, "CAN_cluster2_maf_summary.pdf"))
  plotmafSummary(maf = maf_cluster2, 
                 rmOutlier = TRUE, 
                 addStat = 'median', 
                 dashboard = TRUE)
  dev.off()
  
  return(list(cluster1 = maf_cluster1, cluster2 = maf_cluster2))
}

maf_summaries <- generate_maf_summaries(merged_maf, mutation_cluster_data, OUTPUT_DIR)

# 11.5 Analyze variants per sample
analyze_variants_per_sample <- function(maf_data_list, cluster_data) {
  cat("Analyzing variants per sample...\n")
  
  variants_cluster1 <- maf_data_list$cluster1@variants.per.sample
  variants_cluster2 <- maf_data_list$cluster2@variants.per.sample
  
  # Statistical test
  if (nrow(variants_cluster1) > 0 && nrow(variants_cluster2) > 0) {
    variant_test <- t.test(variants_cluster1$Variants, variants_cluster2$Variants)
    cat(sprintf("Variant count t-test p-value: %.4f\n", variant_test$p.value))
  }
  
  return(list(cluster1 = variants_cluster1, cluster2 = variants_cluster2))
}

variant_analysis <- analyze_variants_per_sample(maf_summaries, mutation_cluster_data)

# 11.6 Calculate Tumor Mutational Burden (TMB)
calculate_tmb <- function(maf_data_list, capture_size = 38) {
  cat("Calculating Tumor Mutational Burden...\n")
  
  tmb_cluster1 <- tmb(maf_data_list$cluster1, captureSize = capture_size, logScale = TRUE)
  tmb_cluster2 <- tmb(maf_data_list$cluster2, captureSize = capture_size, logScale = TRUE)
  
  # Statistical test
  tmb_test <- wilcox.test(tmb_cluster1$total_perMB, tmb_cluster2$total_perMB)
  cat(sprintf("TMB Wilcoxon test p-value: %.4f\n", tmb_test$p.value))
  
  # Prepare TMB data
  tmb_cluster1$cluster <- "1"
  tmb_cluster2$cluster <- "2"
  tmb_combined <- rbind(tmb_cluster1, tmb_cluster2)
  colnames(tmb_combined)[3] <- "TMB"
  
  # Remove outliers (top 5%)
  tmb_filtered <- tmb_combined %>%
    filter(TMB < quantile(TMB, 0.95))
  
  return(tmb_filtered)
}

tmb_results <- calculate_tmb(maf_summaries)

# 11.7 Visualize TMB results
visualize_tmb <- function(tmb_data, output_file) {
  cat("Creating TMB visualization...\n")
  
  tmb_plot <- ggplot(tmb_data, aes(x = cluster, y = TMB, color = cluster)) +
    geom_boxplot(width = 0.5, size = 2, alpha = 0.2, notch = TRUE, notchwidth = 0.5) +
    scale_color_manual(values = CLUSTER_COLORS) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    ylab("TMB score") +
    xlab(NULL) +
    guides(scale = "none") +
    geom_signif(comparisons = list(c("1", "2")),
                map_signif_level = TRUE,
                textsize = 8, color = "black") +
    geom_jitter(width = 0.1, size = 1.5, aes(fill = cluster),
                shape = 21, alpha = 0.6)
  
  ggsave(output_file, tmb_plot, width = 10, height = 10)
  
  return(tmb_plot)
}

tmb_visualization <- visualize_tmb(tmb_results, 
                                   file.path(OUTPUT_DIR, "CAN_TMB_plot.pdf"))

cat("Mutation analysis completed.\n")