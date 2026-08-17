# =========================================================
# WGCNA: (GSE76895)
# Platform: GPL570
# =========================================================

#---------------------------------------------------------
# 1. Setup environment
#---------------------------------------------------------
setwd("E:/NN project 502/WGCNA")

library(GEOquery)
library(WGCNA)
library(limma)
library(tidyverse)

options(stringsAsFactors=FALSE)
allowWGCNAThreads()
#---------------------------------------------------------
# 2. Download GEO dataset
#---------------------------------------------------------
gse <- getGEO("GSE76895", GSEMatrix = TRUE)
gse <- gse[[1]]   # extract ExpressionSet object

#---------------------------------------------------------
# 3. Extract phenotype (sample metadata)
#---------------------------------------------------------
pheno <- pData(gse)

# create group variable
pheno$group <- pheno$`diabetes status (nd (non-diabetic), t2d (type 2 diabetic), igt (impaired glucose tolerance, t3cd (type 3 diabetic)):ch1`

# check sample distribution
table(pheno$group)

# Exclude T3cD samples
pheno <- pheno[!pheno$group %in% c("T3cD"), ]

# convert to factor
pheno$group <- as.factor(pheno$group)

# extract selected sample IDs
selected_samples <- rownames(pheno)


#---------------------------------------------------------
# 4. Extract expression matrix
#---------------------------------------------------------
expr <- exprs(gse)

# keep only selected samples
expr_subset <- expr[, selected_samples]

dim(expr_subset)


#---------------------------------------------------------
# 5. Download platform annotation (GPL570)
#---------------------------------------------------------
gpl <- getGEO("GPL570", AnnotGPL = TRUE)
gpl_table <- Table(gpl)


#---------------------------------------------------------
# 6. Prepare expression dataframe for merging
#---------------------------------------------------------
expr_df <- as.data.frame(expr_subset)
expr_df$ProbeID <- rownames(expr_df)

# select annotation columns
annotation <- gpl_table[, c("ID", "Gene symbol", "Gene title")]
colnames(annotation) <- c("ProbeID", "GeneSymbol", "GeneTitle")


#---------------------------------------------------------
# 7. Merge probe annotation with expression data
#---------------------------------------------------------
expr_annotated <- merge(annotation, expr_df, by = "ProbeID")

# remove probe column
expr <- expr_annotated[, -1]

# rename gene column
colnames(expr)[1] <- "Gene"

# remove GeneTitle column
expr <- expr[, -2]

#---------------------------------------------------------
# 8. Collapse duplicate probes by gene (mean expression)
#---------------------------------------------------------
expr_gene <- expr %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# remove missing gene names
expr_gene <- as.data.frame(expr_gene)
expr_gene <- expr_gene[!is.na(expr_gene$Gene) & expr_gene$Gene != "", ]

str(expr_gene)
#prepare data for WGCNA

rownames(expr_gene) <- expr_gene$Gene
expr_gene <- expr_gene[, -1]
write.csv(expr_gene,"expression_matrix_83 samples.csv")
datExpr<- t(expr_gene)

#detect bad genes

gsg<- goodSamplesGenes(datExpr , verbose = 3)
summary(gsg)
gsg$allOK

#sample clustering(outlier detection)

sampleTree<- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering", sub="", xlab="")

#soft threshold selection
powers<- c(1:20)
sft<- pickSoftThreshold(datExpr, powerVector = powers,
                        verbose = 5,
                        networkType = "signed")

plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit",
     type="n"
    )

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h = 0.8, col = "blue", lty = 2)

#mean connectivity plot

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n")

text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels = powers,
     col = "red")
softpower<- 12

#network construction
adjacency<- adjacency(datExpr , power = softpower,
                      type = "signed")

TOM <- TOMsimilarity(adjacency, TOMType = "signed")  # signed TOM
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity",
     labels=FALSE, hang=0.04)
Modules <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                         deepSplit = 2, pamRespectsDendro = FALSE,
                         minClusterSize = 30)

table(Modules)  # number of genes in each module
ModuleColors <- labels2colors(Modules)  # module number → color
table(ModuleColors)  # number of genes in each color/module
plotDendroAndColors(geneTree, ModuleColors, "Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors (signed network)")


#module eigengene identification

MElist <- moduleEigengenes(datExpr, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

#module merging

ME.dissimilarity <- 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity
METree <- hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree, main = "Clustering of module eigengenes")
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75

merge <- mergeCloseModules(datExpr, ModuleColors, cutHeight = .25)
# The merged module colors, assigning one color to each module
mergedColors <- merge$colors
# Eigengenes of the new merged modules
mergedMEs <- merge$newMEs

#dendrogram with merged & unmerged module colors
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")

#trait data

#order the groups according to progression
pheno$group <- factor(pheno$group, levels = c("ND", "IGT", "T2D"))  # choose ND as reference

#add a numeric column which indicates severity
pheno$Severity <- as.numeric(pheno$group)

#prepare model matrix
traitData <- model.matrix(~ 0 + group,
                          data = pheno)  # remove intercept
#traitData <- model.matrix(~ 0 + group + Severity,data = pheno)
colnames(traitData)

#remove the word 'group' from column name
colnames(traitData) <- gsub("group", "", colnames(traitData))

# Module-trait correlations
# calculate correlations between module eigengenes and traits

#for pearson correlation between ME & trait data
moduleTraitCor <- cor(mergedMEs, traitData, use = "p")  

#calculate P value
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))  # p-values

# View correlation matrix
round(moduleTraitCor, 2)
round(moduleTraitP, 3)

# Plot labeled heatmap
sizeGrWindow(10,6)
par(mar = c(10, 15, 3, 3))
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitP, 1), ")", sep = "")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traitData),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-trait relationships (signed network)")

save.image("wgcna_updated.RData")

