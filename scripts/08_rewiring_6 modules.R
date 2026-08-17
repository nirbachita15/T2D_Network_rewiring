

library(WGCNA)
library(tidyverse)

#=========================================================
# 1. Selected genes (combined modules)
#=========================================================

selected_modules <- c(
  "brown",
  "darkgreen",
  "lightcyan",
  "midnightblue",
  "orange",
  "salmon"
)

genes_all <- colnames(datExpr)[mergedColors %in% selected_modules]
genes_all <- unique(genes_all)

expr_all <- datExpr[, genes_all]

#=========================================================
# 2. Split groups
#=========================================================

ND  <- expr_all[pheno$group == "ND", ]
IGT <- expr_all[pheno$group == "IGT", ]
T2D <- expr_all[pheno$group == "T2D", ]

#=========================================================
# 3. Correlation matrices
#=========================================================

cor_ND  <- cor(ND,  use = "pairwise.complete.obs")
cor_IGT <- cor(IGT, use = "pairwise.complete.obs")
cor_T2D <- cor(T2D, use = "pairwise.complete.obs")

#=========================================================
# 4. FAST edge rewiring function
#=========================================================

get_rewired_edges <- function(cor1, cor2, g1, g2, threshold = 0.7) {
  
  idx <- upper.tri(cor1)
  
  df <- data.frame(
    Gene1 = rownames(cor1)[row(cor1)[idx]],
    Gene2 = colnames(cor1)[col(cor1)[idx]],
    Cor1  = cor1[idx],
    Cor2  = cor2[idx]
  )
  
  df$Diff <- abs(df$Cor1 - df$Cor2)
  
  # keep only strong rewiring edges
  df <- df %>% filter(Diff > threshold)
  
  colnames(df)[3] <- paste0(g1, "_cor")
  colnames(df)[4] <- paste0(g2, "_cor")
  
  df <- df %>% arrange(desc(Diff))
  
  return(df)
}

#=========================================================
# 5. Run comparisons
#=========================================================

ND_vs_IGT  <- get_rewired_edges(cor_ND, cor_IGT, "ND", "IGT", threshold = 0.7)
ND_vs_T2D  <- get_rewired_edges(cor_ND, cor_T2D, "ND", "T2D", threshold = 0.7)
IGT_vs_T2D <- get_rewired_edges(cor_IGT, cor_T2D, "IGT", "T2D", threshold = 0.7)

#=========================================================
# 6. Save outputs
#=========================================================

dir.create("rewired_edges", showWarnings = FALSE)

write.csv(ND_vs_IGT,
          "rewired_edges/ND_vs_IGT_edges.csv",
          row.names = FALSE)

write.csv(ND_vs_T2D,
          "rewired_edges/ND_vs_T2D_edges.csv",
          row.names = FALSE)

write.csv(IGT_vs_T2D,
          "rewired_edges/IGT_vs_T2D_edges.csv",
          row.names = FALSE)

#=========================================================
# 7. Check results
#=========================================================

head(ND_vs_IGT)
head(ND_vs_T2D)
head(IGT_vs_T2D)


getwd()

