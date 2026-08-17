# =========================================================
# Differential Expression Analysis: IGT vs ND (GSE76895)
# Platform: GPL570
# =========================================================

#---------------------------------------------------------
# 1. Setup environment
#---------------------------------------------------------
setwd("E:/NN project 502/DEG/dhrubo_DEG")

library(GEOquery)
library(tidyverse)
library(limma)

#---------------------------------------------------------
# 2. Download GEO dataset
#---------------------------------------------------------
gse <- getGEO("GSE76895", GSEMatrix = TRUE)
gse <- gse[[1]]   # extract ExpressionSet object

experimentData(gse)
#---------------------------------------------------------
# 3. Extract phenotype (sample metadata)
#---------------------------------------------------------
pheno <- pData(gse)

write.csv(pheno, "GSE76895_pheno_data.csv")
# create group variable
pheno$group <- pheno$`diabetes status (nd (non-diabetic), t2d (type 2 diabetic), igt (impaired glucose tolerance, t3cd (type 3 diabetic)):ch1`

write.csv(pheno, "pheno_with_group_column.csv")
# check sample distribution
table(pheno$group)

# keep only ND and IGT samples
pheno <- pheno[!pheno$group %in% c("T2D", "T3cD"), ]

write.csv(pheno, "pheno_ND_IGT_only.csv")
# convert to factor
pheno$group <- as.factor(pheno$group)

# extract selected sample IDs
selected_samples <- rownames(pheno)


#---------------------------------------------------------
# 4. Extract expression matrix
#---------------------------------------------------------
expr <- exprs(gse)

write.csv(expr, "raw_expression_matrix.csv")
# keep only selected samples
expr_subset <- expr[, selected_samples]

dim(expr_subset)

# save raw subset matrix
write.csv(expr_subset, "GSE76895_expression_subset.csv", row.names = TRUE)


#---------------------------------------------------------
# 5. Download platform annotation (GPL570)
#---------------------------------------------------------
gpl <- getGEO("GPL570", AnnotGPL = TRUE)
gpl_table <- Table(gpl)

write.csv(gpl_table, "GPL570_annotation_full.csv")

#---------------------------------------------------------
# 6. Prepare expression dataframe for merging
#---------------------------------------------------------
expr_df <- as.data.frame(expr_subset)
expr_df$ProbeID <- rownames(expr_df)

write.csv(expr_df, "expression_with_probeID.csv")
# select annotation columns
annotation <- gpl_table[, c("ID", "Gene symbol", "Gene title")]
colnames(annotation) <- c("ProbeID", "GeneSymbol", "GeneTitle")

write.csv(annotation, "probe_gene_annotation.csv")

#---------------------------------------------------------
# 7. Merge probe annotation with expression data
#---------------------------------------------------------
expr_annotated <- merge(annotation, expr_df, by = "ProbeID")
write.csv(expr_annotated, "expression_with_annotation.csv")

# remove probe column
expr <- expr_annotated[, -1]

# rename gene column
colnames(expr)[1] <- "Gene"

# remove GeneTitle column
expr <- expr[, -2]

write.csv(expr, "expression_gene_symbol_table.csv")
#---------------------------------------------------------
# 8. Collapse duplicate probes by gene (mean expression)
#---------------------------------------------------------
expr_gene <- expr %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# remove missing gene names
expr_gene <- as.data.frame(expr_gene)
expr_gene <- expr_gene[!is.na(expr_gene$Gene) & expr_gene$Gene != "", ]

write.csv(expr_gene, "gene_level_expression_clean.csv")
#---------------------------------------------------------
# 9. Prepare expression matrix for limma
#---------------------------------------------------------
rownames(expr_gene) <- expr_gene$Gene
expr_gene <- expr_gene[, -1]

expr_matrix <- as.matrix(expr_gene)

write.csv(expr_matrix, "expression_matrix_for_limma.csv")
#---------------------------------------------------------
# 10. Prepare sample group vector
#---------------------------------------------------------
group <- pheno$group
names(group) <- rownames(pheno)

group <- factor(group)
table(group)


#---------------------------------------------------------
# 11. Create design matrix
#---------------------------------------------------------
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

design


#---------------------------------------------------------
# 12. Fit linear model using limma
#---------------------------------------------------------
fit <- lmFit(expr_matrix, design)


#---------------------------------------------------------
# 13. Define contrast (IGT vs ND)
#---------------------------------------------------------
contrast_matrix <- makeContrasts(
  IGT_vs_ND = IGT - ND,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


#---------------------------------------------------------
# 14. Extract differential expression results
#---------------------------------------------------------
deg_results <- topTable(
  fit2,
  coef = "IGT_vs_ND",
  number = Inf,
  adjust.method = "BH"
)

View(deg_results)

write.csv(deg_results, "IGT_vs_ND_DEG_results.csv")
# number of genes with raw p < 0.05
sum(deg_results$P.Value < 0.05)


upregulated <- deg_results[deg_results$logFC > 1 & deg_results$P.Value < 0.05, ]

downregulated <- deg_results[deg_results$logFC < -1 & deg_results$P.Value < 0.05, ]

write.csv(upregulated, "upregulated_genes.csv")
write.csv(downregulated, "downregulated_genes.csv")







