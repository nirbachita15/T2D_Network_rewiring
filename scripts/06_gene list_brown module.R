# 1. datExpr er matrix theke column name (Gene Symbols) evong updated merged colors alada kora
all_genes <- colnames(datExpr)

# 2. Sudhu matro 'brown' module er gene gulo ke filter kore kete neya
brown_module_genes <- all_genes[mergedColors == "brown"]

# 3. Dataframe e convert kora jate dekhte o transfer korte subidha hoy
brown_genes_df <- data.frame(Gene = brown_module_genes)

# 4. Apnar set-working-directory (setwd) te excel/csv file hisebe seve kora
write.csv(brown_genes_df, file = "WGCNA_Merged_Brown_Module_Genes.csv", row.names = FALSE)

# 5. Check kora total koyti gene eii brown module e asche
print(paste("Total genes in Merged Brown Module:", length(brown_module_genes)))
