#GO and kegg
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# 1. Take all the gene names from the column of the expression matrix
all_genes <- colnames(datExpr)

# ২. filter only brown module gene 
brown_genes <- all_genes[mergedColors == "brown"]

# check the brown module gene 
length(brown_genes)

head(brown_genes)
gene_ids <- bitr(brown_genes, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)
cat("Mapped:", nrow(gene_ids), "out of", length(brown_genes), "\n")



head(gene_ids)


# Run GO Enrichment 
go_BP <- enrichGO(    gene          = gene_ids$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",        # BP
                      pAdjustMethod = "BH",         # Benjamini-Hochberg p-value adjustment
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)         
go_CC <- enrichGO(    gene          = gene_ids$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          ont           = "CC",        # CC
                          pAdjustMethod = "BH",         # Benjamini-Hochberg p-value adjustment
                          pvalueCutoff  = 0.1,
                          qvalueCutoff  = 0.1,
                          readable      = TRUE)   
go_MF <- enrichGO(    gene          = gene_ids$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "MF",        #  MF 
                      pAdjustMethod = "BH",         # Benjamini-Hochberg p-value adjustment
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.1,
                      readable      = TRUE)         
# View results as dataframes
go_BP_df <- as.data.frame(go_BP)
go_CC_df <- as.data.frame(go_CC)
go_MF_df <- as.data.frame(go_MF)


# Save results
write.csv(go_BP_df, "Brown_module_GO_BP.csv", row.names = FALSE)
write.csv(go_CC_df, "Brown_module_GO_CC.csv", row.names = FALSE)
write.csv(go_MF_df, "Brown_module_GO_MF.csv", row.names = FALSE)

#barplot

barplot(go_BP, showCategory = 10, title = "GO Biological Process - Brown Module")
barplot(go_CC, showCategory = 10, title = "GO Cellular Component - Brown Module")
barplot(go_MF, showCategory = 10, title = "GO Molecular Function - Brown Module")

p_bp <- barplot(go_BP, showCategory = 10, title = "GO Biological Process - Brown Module")
p_cc <- barplot(go_CC_new, showCategory = 10, title = "GO Cellular Component - Brown Module")
p_mf <- barplot(go_MF, showCategory = 10, title = "GO Molecular Function - Brown Module")
# save the plot
ggsave("GO_BP_brown_module.png", plot = p_bp, width = 8, height = 6, dpi = 300)
ggsave("GO_CC_brown_module1.png", plot = p_cc, width = 8, height = 6, dpi = 300)
ggsave("GO_MF_brown_module.png", plot = p_mf, width = 8, height = 6, dpi = 300)


# KEGG Enrichment (for Human , organism  = "hsa")
kegg_enrich <- enrichKEGG(gene          = gene_ids$ENTREZID,
                          organism      = "hsa",     
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.1,
                          qvalueCutoff  = 0.1)

# Convert the result into Gene symbol 
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# save the result as  excel file (.csv)
write.csv(as.data.frame(kegg_enrich), "Brown_Module_KEGG_Results.csv", row.names = FALSE)

# barplot 
barplot(kegg_enrich,showCategory = 10,title = "KEGG Pathway Enrichment for Brown module")
p_kegg <- barplot(kegg_enrich,showCategory = 10, title = "KEGG Pathway Enrichment for Brown module")
ggsave("KEGG_Brown_Module.png", plot = p_kegg, width = 8, height = 6, dpi = 300)


#---------for supplimentary file----------
# total significant term
num_significant_terms_BP <- nrow(as.data.frame(go_BP))
num_significant_terms_CC <- nrow(as.data.frame(go_CC_new))
num_significant_terms_MF <- nrow(as.data.frame(go_MF))

# barplot- showing all the significant term
barplot(go_BP, 
        showCategory = num_significant_terms_BP, 
        title = "GO Biological Process - Brown Module")
barplot(go_CC, 
        showCategory = num_significant_terms_CC, 
        title = "GO Cellular Component - Brown Module")
barplot(go_MF, 
        showCategory = num_significant_terms_MF, 
        title = "GO Molecular Function - Brown Module")
