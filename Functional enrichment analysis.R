# ============================================================
# Functional Enrichment Analysis
# ============================================================

cat("Performing functional enrichment analysis...\n")

#Prepare gene list for enrichment analysis
cat("Preparing gene list for enrichment analysis...\n")
gene_symbols <- selected_genes

# Convert gene symbols to Entrez IDs
convert_gene_ids <- function(gene_symbols) {
  id_conversion <- bitr(gene_symbols,
                        fromType = "SYMBOL",
                        toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                        OrgDb = "org.Hs.eg.db")
  
  id_conversion <- na.omit(id_conversion)
  return(id_conversion)
}

gene_id_list <- convert_gene_ids(gene_symbols)

#KEGG pathway enrichment analysis
perform_kegg_enrichment <- function(gene_ids, output_dir) {
  cat("Performing KEGG pathway enrichment analysis...\n")
  
  kegg_results <- enrichKEGG(gene = gene_ids$ENTREZID,
                             organism = "hsa",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)
  
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    kegg_dataframe <- data.frame(kegg_results)
    
    # Create visualizations
    kegg_dotplot <- dotplot(kegg_results, orderBy = "x", showCategory = 20)
    kegg_barplot <- barplot(kegg_results, orderBy = "x", showCategory = 20)
    
    # Save plots
    ggsave(file.path(output_dir, "CAN_kegg_dotplot.pdf"), kegg_dotplot, width = 10, height = 10)
    ggsave(file.path(output_dir, "CAN_kegg_barplot.pdf"), kegg_barplot, width = 10, height = 10)
    
    return(kegg_results)
  } else {
    cat("No significant KEGG pathways found.\n")
    return(NULL)
  }
}

kegg_enrichment <- perform_kegg_enrichment(gene_id_list, OUTPUT_DIR)

cat("Functional enrichment analysis completed.\n")