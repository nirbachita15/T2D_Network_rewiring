# install package

Sys.which("make")

install.packages(c(
  "devtools",
  "remotes",
  "stringi",
  "stringr"
))

install.packages("BiocManager")
BiocManager::install(c(
  "qvalue",
  "impute",
  "GO.db",
  "AnnotationDbi"
))

remotes::install_github(
  "andymckenzie/DGCA",
  dependencies = TRUE,
  upgrade = "never"
)
# =========================================================
# MODULE-BASED NETWORK REWIRING ANALYSIS
# =========================================================

#---------------------------------------------------------
# Load packages
#---------------------------------------------------------

library(DGCA)
library(WGCNA)
library(dplyr)
library(tidyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)

dim(datExpr_ND)
dim(datExpr_IGT)
dim(datExpr_T2D)

#check gene order----
all(colnames(datExpr_ND)==colnames(datExpr_T2D))

all(colnames(datExpr_IGT)==colnames(datExpr_T2D))

#select module
selectedmodules<- "brown"

#extract gene from module
selectedgenes<- modulecolors%in%selectedmodules
genesuse<- colnames(datExpr_ND)[selectedgenes]
length(genesuse)


expr_nd<- datExpr_ND[,genesuse]
expr_t2d<- datExpr_T2D[,genesuse]
expr_igt<- datExpr_IGT[,genesuse]

#transpose matrix----
expr_nd<- t(expr_nd)
expr_igt<- t(expr_igt)
expr_t2d<- t(expr_t2d)

#check dimension of new created expression matrix----
dim(expr_nd)
dim(expr_igt)
dim(expr_t2d)

#combine expression matrix----
expr_combined1<- cbind(expr_nd,expr_t2d)

#create group vector----
groupVec<- c(rep("ND",ncol(expr_nd)),
             rep("T2D",ncol(expr_t2d)))

design<- makeDesign(groupVec)

ddcResults1<- ddcorAll(
  inputMat = expr_combined1,
  design = design,
  compare = c("ND","T2D"),
  adjust = "perm",
  nPerms = 100,
  corrType = "pearson",
  heatmapPlot = FALSE
)
head(ddcResults1)

ddcResults1<- ddcResults1%>%
  mutate(dcor = T2D_cor - ND_cor)

sigEdges<- ddcResults1%>%
  filter(pValDiff<0.05,
         abs(dcor)>0.8)

nrow(sigEdges)

nodes <- sigEdges %>%
  
  pivot_longer(
    cols = c(Gene1, Gene2),
    names_to = "end",
    values_to = "gene"
  ) %>%
  
  group_by(gene) %>%
  
  summarise(
    
    rewiring_score = sum(abs(dcor)),
    
    rewired_degree = n(),
    
    .groups = "drop"
  )


top_rewired <- nodes %>%
  
  arrange(desc(rewiring_score))

head(top_rewired, 20)

write.csv(
  sigEdges,
  "T2D_vs_ND_rewired_edges.csv",
  row.names = FALSE
)

write.csv(
  top_rewired,
  "T2D_vs_ND_rewired_nodes.csv",
  row.names = FALSE
)

# =========================================================
# Clean file for cytoscape
# =========================================================
library(dplyr)
library(stringr)
library(tidyr)

sigEdges_clean <- sigEdges %>%
  mutate(
    Gene1 = str_split(Gene1, "///") %>% sapply("[", 1),
    Gene2 = str_split(Gene2, "///") %>% sapply("[", 1)
  )


nodes_clean <- sigEdges_clean %>%
  pivot_longer(
    cols = c(Gene1, Gene2),
    names_to = "end",
    values_to = "ID"  
  ) %>%
  group_by(ID) %>%
  summarise(
    rewiring_score = sum(abs(dcor)),
    rewired_degree = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(rewiring_score))


write.csv(sigEdges_clean, "Edges_Clean.csv", row.names = FALSE)
write.csv(nodes_clean, "Nodes_Clean.csv", row.names = FALSE)



library(dplyr)
library(stringr)
library(tidyr)

library(dplyr)
library(stringr)
library(tidyr)

edges_clean <- sigEdges %>%
  mutate(
    Gene1 = str_split(Gene1, "///") %>% sapply("[", 1) %>% str_trim(),
    Gene2 = str_split(Gene2, "///") %>% sapply("[", 1) %>% str_trim()
  )


scores_g1 <- edges_clean %>%
  pivot_longer(cols = c(Gene1, Gene2), names_to = "type", values_to = "Gene") %>%
  group_by(Gene) %>%
  summarise(
    Gene1_Score = sum(abs(dcor)),   
    Gene1_Degree = n(),
    .groups = "drop"
  )


scores_g2 <- edges_clean %>%
  pivot_longer(cols = c(Gene1, Gene2), names_to = "type", values_to = "Gene") %>%
  group_by(Gene) %>%
  summarise(
    Gene2_Score = sum(abs(dcor)),    
    Gene2_Degree = n(),
    .groups = "drop"
  )


edges_with_node_data <- edges_clean %>%
  left_join(scores_g1, by = c("Gene1" = "Gene")) %>%
  left_join(scores_g2, by = c("Gene2" = "Gene"))

write.csv(edges_with_node_data, "All_In_One_Network.csv", row.names = FALSE)
