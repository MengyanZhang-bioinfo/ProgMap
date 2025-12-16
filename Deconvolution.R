# ============================================================
# Deconvolution Analysis
# ============================================================

cat("Performing deconvolution analysis...\n")

# Run CIBERSORT for original stages
cat("Running CIBERSORT for original stages...\n")
setwd(STAGE1_DIR)
source("cibersort.R")
CIBERSORT("LM4_mean.txt", "CAN_tpm_stage1.txt", 
          perm = PERMUTATIONS, QN = TRUE, 
          "stage1_CIBERSORT-Results_LM4_mean.txt")

# Stage 2 analysis
setwd(STAGE2_DIR)
source("cibersort.R")
CIBERSORT("LM4_mean.txt", "CAN_tpm_stage2.txt", 
          perm = PERMUTATIONS, QN = TRUE, 
          "stage2_CIBERSORT-Results_LM4_mean.txt")

# 5.2 Load and process results
stage1_ciber <- read.table(file.path(STAGE1_DIR, "stage1_CIBERSORT-Results_LM4_mean.txt"), 
                           header = TRUE, sep = "\t")
stage2_ciber <- read.table(file.path(STAGE2_DIR, "stage2_CIBERSORT-Results_LM4_mean.txt"), 
                           header = TRUE, sep = "\t")

# Combine results
combined_ciber <- cbind(t(stage1_ciber[, 1:(ncol(stage1_ciber) - 3)]), 
                        t(stage2_ciber[, 1:(ncol(stage1_ciber) - 3)]))

colnames(combined_ciber) <- combined_ciber[1, ]
combined_ciber <- combined_ciber[-1, ]
combined_ciber <- as.data.frame(combined_ciber)
combined_ciber[, 1:ncol(combined_ciber)] <- apply(combined_ciber[, 1:ncol(combined_ciber)], 2, as.numeric)

# Statistical analysis
p_values <- numeric(nrow(combined_ciber))
for (j in 1:nrow(combined_ciber)) {
  test_result <- wilcox.test(as.numeric(combined_ciber[j, 1:nrow(stage1_ciber)]), 
                             as.numeric(combined_ciber[j, (nrow(stage1_ciber) + 1):ncol(combined_ciber)]), 
                             paired = FALSE, var.equal = TRUE)
  p_values[j] <- test_result$p.value
}

# 5.5 Visualization
plot_data <- data.frame()
cell_types <- rep(c("stage1", "stage2"), c(nrow(stage1_ciber), nrow(stage2_ciber)))

for (i in 1:nrow(combined_ciber)) {
  scores <- as.numeric(combined_ciber[i, ])
  cell_type <- rep(rownames(combined_ciber)[i], ncol(combined_ciber))
  stage <- cell_types
  plot_data <- rbind(plot_data, data.frame(scores, cell_type, stage))
}

plot_data$scores <- as.numeric(plot_data$scores)
plot_data$stage <- as.factor(plot_data$stage)

# Create visualization
stage_deconvolution_plot <- ggplot(plot_data, aes(y = scores, x = factor(cell_type), fill = stage)) +
  geom_boxplot(outlier.size = 1) +
  theme_classic() +
  xlab("") + ylab("") +
  scale_fill_manual(values = STAGE_COLORS) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")

ggsave(file.path(OUTPUT_DIR, "CAN_stage_deconvolution_plot.pdf"), 
       stage_deconvolution_plot, width = 10, height = 10)

cat("Deconvolution analysis completed.\n")