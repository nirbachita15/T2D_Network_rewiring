#hub gene identification----
#use merged colors as module colors----
modulecolors<- mergedColors
#recalculate eigengene using merged modules----
MEs0 <- moduleEigengenes(datExpr, modulecolors)$eigengenes

MEs <- orderMEs(MEs0)

#extract module names----
modNames <- substring(names(MEs), 3)

#choose stage----
targetTrait <- "T2D"

#=========================================================
# CREATE BINARY TRAIT----
#=========================================================

trait <- as.data.frame(traitData[, targetTrait])

names(trait) <- targetTrait

head(trait)

#=========================================================
# CALCULATE MODULE MEMBERSHIP (MM)----
#=========================================================

geneModuleMembership <- as.data.frame(
  cor(datExpr, MEs, use = "p")
)

MMPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(geneModuleMembership),
    nrow(datExpr)
  )
)

names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

#=========================================================
# CALCULATE GENE SIGNIFICANCE (GS)----
#=========================================================

geneTraitSignificance <- as.data.frame(
  cor(datExpr, trait, use = "p")
)

GSPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(geneTraitSignificance),
    nrow(datExpr)
  )
)

names(geneTraitSignificance) <- "GS.T2D"
names(GSPvalue) <- "p.GS.T2D"

#for T2D hub genes-----
# SELECT MODULES OF INTEREST----
#=========================================================

selectedModules <- c(
  "darkgreen",
  "salmon",
  "brown",
  "midnightblue",
  "orange", "lightcyan"
)

# EMPTY LIST TO STORE RESULTS----
#=========================================================

allHubGenes <- list()

#=========================================================
# LOOP THROUGH MODULES----
#=========================================================

for(module in selectedModules){
  
  print(paste("Processing module:", module))
  
  # find module column
  column <- match(module, modNames)
  
  # identify genes in module
  moduleGenes <- modulecolors == module
  #=======================================================
  # CREATE HUB GENE TABLE
  #=======================================================
  
  hubGenes <- data.frame(
    
    Gene = colnames(datExpr)[moduleGenes],
    
    MM = geneModuleMembership[moduleGenes, column],
    
    MM_pvalue = MMPvalue[moduleGenes, column],
    
    GS = geneTraitSignificance[moduleGenes, 1],
    
    GS_pvalue = GSPvalue[moduleGenes, 1],
    
    Module = module
  )
  
  #=======================================================
  # FILTER HUB GENES
  #=======================================================
  
  hubGenes_filtered <- hubGenes[
    
    abs(hubGenes$MM) > 0.7 &
      abs(hubGenes$GS) > 0.2,
    
  ]
  
  # sort by MM
  hubGenes_filtered <- hubGenes_filtered[
    order(-abs(hubGenes_filtered$MM)),
  ]
  
  #=======================================================
  # SAVE EACH MODULE RESULT
  #=======================================================
  
  write.csv(
    hubGenes_filtered,
    
    paste0("HubGenes_", module, "_T2D.csv"),
    
    row.names = FALSE
  )
  
  #=======================================================
  # STORE IN LIST
  #=======================================================
  
  allHubGenes[[module]] <- hubGenes_filtered
  
  #=======================================================
  # OPTIONAL SCATTERPLOT
  #=======================================================
  
  verboseScatterplot(
    
    abs(geneModuleMembership[moduleGenes, column]),
    
    abs(geneTraitSignificance[moduleGenes, 1]),
    
    xlab = paste("Module Membership in", module),
    
    ylab = "Gene significance for T2D",
    
    main = paste(module,
                 "Module Membership vs Gene Significance"),
    
    col = module
  )
}

#for IGT hub genes----

# SELECT MODULES OF INTEREST
#=========================================================

selectedModules <- c(
  "midnightblue",
  "darkgreen",
  "brown","orange"
)

#=========================================================
# CREATE IGT TRAIT
#=========================================================

trait <- as.data.frame(traitData[, "IGT"])

names(trait) <- "IGT"

#=========================================================
# CALCULATE GENE SIGNIFICANCE FOR IGT
#=========================================================

geneTraitSignificance <- as.data.frame(
  
  cor(datExpr,
      trait,
      use = "p")
)

GSPvalue <- as.data.frame(
  
  corPvalueStudent(
    
    as.matrix(geneTraitSignificance),
    
    nrow(datExpr)
  )
)

names(geneTraitSignificance) <- "GS.IGT"

names(GSPvalue) <- "p.GS.IGT"

#=========================================================
# EMPTY LIST TO STORE RESULTS
#=========================================================

allHubGenes <- list()

#=========================================================
# LOOP THROUGH MODULES
#=========================================================

for(module in selectedModules){
  
  print(paste("Processing module:", module))
  
  #=======================================================
  # FIND MODULE COLUMN
  #=======================================================
  
  column <- match(module, modNames)
  
  #=======================================================
  # IDENTIFY GENES IN MODULE
  #=======================================================
  
  moduleGenes <- modulecolors == module
  
  #=======================================================
  # CREATE HUB GENE TABLE
  #=======================================================
  
  hubGenes <- data.frame(
    
    Gene = colnames(datExpr)[moduleGenes],
    
    MM = geneModuleMembership[moduleGenes, column],
    
    MM_pvalue = MMPvalue[moduleGenes, column],
    
    GS = geneTraitSignificance[moduleGenes, 1],
    
    GS_pvalue = GSPvalue[moduleGenes, 1],
    
    Module = module
  )
  
  #=======================================================
  # FILTER HUB GENES
  #=======================================================
  
  hubGenes_filtered <- hubGenes[
    
    abs(hubGenes$MM) > 0.7 &
      abs(hubGenes$GS) > 0.2,
    
  ]
  
  #=======================================================
  # SORT BY MM
  #=======================================================
  
  hubGenes_filtered <- hubGenes_filtered[
    
    order(-abs(hubGenes_filtered$MM)),
  ]
  
  #=======================================================
  # VIEW TOP HUB GENES
  #=======================================================
  
  print(head(hubGenes_filtered))
  
  #=======================================================
  # SAVE EACH MODULE RESULT
  #=======================================================
  
  write.csv(
    
    hubGenes_filtered,
    
    paste0("HubGenes_", module, "_IGT_Nirbachita.csv"),
    
    row.names = FALSE
  )
  
  #=======================================================
  # STORE IN LIST
  #=======================================================
  
  allHubGenes[[module]] <- hubGenes_filtered
  
  #=======================================================
  # OPTIONAL SCATTERPLOT
  #=======================================================
  
  verboseScatterplot(
    
    abs(geneModuleMembership[moduleGenes, column]),
    
    abs(geneTraitSignificance[moduleGenes, 1]),
    
    xlab = paste("Module Membership in", module),
    
    ylab = "Gene significance for IGT",
    
    main = paste(module,
                 "Module Membership vs Gene Significance"),
    
    col = module
  )
}

#for ND hub genes----
# SELECT MODULES OF INTEREST
#=========================================================

selectedModules <- c(
  "salmon",
  "brown", "lightcyan"
)

#=========================================================
# CREATE ND TRAIT
#=========================================================

trait <- as.data.frame(traitData[, "ND"])

names(trait) <- "ND"

#=========================================================
# CALCULATE GENE SIGNIFICANCE FOR IGT
#=========================================================

geneTraitSignificance <- as.data.frame(
  
  cor(datExpr,
      trait,
      use = "p")
)

GSPvalue <- as.data.frame(
  
  corPvalueStudent(
    
    as.matrix(geneTraitSignificance),
    
    nrow(datExpr)
  )
)

names(geneTraitSignificance) <- "GS.IGT"

names(GSPvalue) <- "p.GS.IGT"

#=========================================================
# EMPTY LIST TO STORE RESULTS
#=========================================================

allHubGenes <- list()

#=========================================================
# LOOP THROUGH MODULES
#=========================================================

for(module in selectedModules){
  
  print(paste("Processing module:", module))
  
  #=======================================================
  # FIND MODULE COLUMN
  #=======================================================
  
  column <- match(module, modNames)
  
  #=======================================================
  # IDENTIFY GENES IN MODULE
  #=======================================================
  
  moduleGenes <- modulecolors == module
  
  #=======================================================
  # CREATE HUB GENE TABLE
  #=======================================================
  
  hubGenes <- data.frame(
    
    Gene = colnames(datExpr)[moduleGenes],
    
    MM = geneModuleMembership[moduleGenes, column],
    
    MM_pvalue = MMPvalue[moduleGenes, column],
    
    GS = geneTraitSignificance[moduleGenes, 1],
    
    GS_pvalue = GSPvalue[moduleGenes, 1],
    
    Module = module
  )
  
  #=======================================================
  # FILTER HUB GENES
  #=======================================================
  
  hubGenes_filtered <- hubGenes[
    
    abs(hubGenes$MM) > 0.7 &
      abs(hubGenes$GS) > 0.2,
    
  ]
  
  #=======================================================
  # SORT BY MM
  #=======================================================
  
  hubGenes_filtered <- hubGenes_filtered[
    
    order(-abs(hubGenes_filtered$MM)),
  ]
  
  #=======================================================
  # VIEW TOP HUB GENES
  #=======================================================
  
  print(head(hubGenes_filtered))
  
  #=======================================================
  # SAVE EACH MODULE RESULT
  #=======================================================
  
  write.csv(
    
    hubGenes_filtered,
    
    paste0("HubGenes_", module, "_ND_Nirbachita.csv"),
    
    row.names = FALSE
  )
  
  #=======================================================
  # STORE IN LIST
  #=======================================================
  
  allHubGenes[[module]] <- hubGenes_filtered
  
  #=======================================================
  # OPTIONAL SCATTERPLOT
  #=======================================================
  
  verboseScatterplot(
    
    abs(geneModuleMembership[moduleGenes, column]),
    
    abs(geneTraitSignificance[moduleGenes, 1]),
    
    xlab = paste("Module Membership in", module),
    
    ylab = "Gene significance for ND",
    
    main = paste(module,
                 "Module Membership vs Gene Significance"),
    
    col = module
  )
}


#=========================================================
# EIGENGENE NETWORK ANALYSIS----
#=========================================================

#---------------------------------------------------------
# CREATE TRAIT DATAFRAME----
#---------------------------------------------------------

Traits <- data.frame(
  
  ND = traitData[, "ND"],
  
  IGT = traitData[, "IGT"],
  
  T2D = traitData[, "T2D"]
)

head(Traits)

#=========================================================
# ADD TRAITS TO MODULE EIGENGENES----
#=========================================================

MET <- orderMEs(
  
  cbind(
    
    mergedMEs,
    
    Traits
  )
)

#=========================================================
# PLOT EIGENGENE NETWORK----
# DENDROGRAM + HEATMAP----
#=========================================================

par(cex = 0.9)

plotEigengeneNetworks(
  
  MET,
  
  "Eigengene Network Visualization",
  
  marDendro = c(0,4,1,2),
  
  marHeatmap = c(7,4,1,2),
  
  cex.lab = 0.8,
  
  xLabelsAngle = 90
)

#=========================================================
# PLOT ONLY DENDROGRAM----
#=========================================================

par(cex = 1.0)

plotEigengeneNetworks(
  
  MET,
  
  "Eigengene dendrogram",
  
  marDendro = c(0,4,2,0),
  
  plotHeatmaps = FALSE
)

#=========================================================
# PLOT ONLY HEATMAP----
#=========================================================

par(cex = 1.0,
    mar = c(1,1,1,1))

plotEigengeneNetworks(
  
  MET,
  
  "Eigengene adjacency heatmap",
  
  marHeatmap = c(5,5,2,2),
  
  plotDendrograms = FALSE,
  
  xLabelsAngle = 90
)


