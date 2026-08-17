
# Download GSE76895
options(timeout = 6000)
gse <- getGEO("GSE76895", GSEMatrix = TRUE)
save.image(file="gse")
gse <- gse[[1]] 

#---------------------------------------------------------
# 3. Extract phenotype (sample metadata)
#---------------------------------------------------------
pheno <- pData(gse)

# Extract the diabetes status column
pheno$group <- pheno[,"diabetes status (nd (non-diabetic), t2d (type 2 diabetic), igt (impaired glucose tolerance, t3cd (type 3 diabetic)):ch1"]

# keep only ND
pheno <- pheno[!pheno$group %in% c("IGT", "T2D", "T3cD"), ]

# convert to factor
pheno$group <- as.factor(pheno$group)

# extract selected sample IDs
selected_samples <- rownames(pheno)

# Extract expression matrix
#---------------------------------------------------------
expr <- exprs(gse)

# keep only selected samples
expr_subset <- expr[, selected_samples]

dim(expr_subset)


# 5. Download platform annotation (GPL570)
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


# 9. Filter out any rows where the Gene name is missing or empty
expr_gene <- as.data.frame(expr_gene)
expr_gene <- expr_gene[!is.na(expr_gene$Gene) & expr_gene$Gene != "", ]

# 10. Set the Gene Symbols as the row names
rownames(expr_gene) <- expr_gene$Gene
expr_gene <- expr_gene[, -1]
datExpr0 <- t(expr_gene)

write.csv(datExpr0, file = "datExpr0_transposed.csv")


# 2. Check for genes and samples with too many missing values
# WGCNA has a specific function for this cleaning step
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
summary(gsg)
gsg$allOK

# Choose a set of soft-thresholding powers

powers<- c(1:20)
spt<- pickSoftThreshold(datExpr0, powerVector = powers,
                        verbose = 5,
                        networkType = "signed")

spt
#plot this data frame to better visualize what β value you should choose.

pdf("Scale_Independence_&_Mean_Connectivity_Plot_ND.pdf", width = 8, height = 6)
#Plotting the results 
par(mfrow = c(1, 2))  
par(mar = c(4.5, 4.5, 2.5, 1.5))

#Index the scale-free topology adjustment as a function of the power soft thresholding. 

plot(spt$fitIndices[,1],-sign(spt$fitIndices[,3])*spt$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", 
     main = paste("Scale independence for ND")) 
text(spt$fitIndices[,1], -sign(spt$fitIndices[,3])*spt$fitIndices[,2], 
     labels=powers,cex=0.9,col="red") 

#This line corresponds to use a cut-off R2 of h 
abline(h=0.8,col="red") 

#Connectivity means as a function of soft power thresholding 
plot(spt$fitIndices[,1], spt$fitIndices[,5], 
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", 
     main = paste("Mean connectivity for ND")) 
text(spt$fitIndices[,1], spt$fitIndices[,5], labels=powers, 
     cex=0.9,col="red") 

#This line corresponds to use a cut-off R2 of h 
abline(h=0.8,col="red")

dev.off()
graphics.off()
#Calling the Adjacency Function
softPower <- 13
adjacency <- adjacency(datExpr0, power = softPower, type = "signed")

netwk_ND <- blockwiseModules(
  datExpr0, 
  power = 13,                # The optimal power we found!
  corType = "bicor",         # If using biweight midcorrelation
  maxBlockSize = 10000,
  # Force WGCNA to use its own internal bicor options:
  networkType = "signed",    
  TOMType = "signed",
  verbose = 3
)

mergedColors = labels2colors(netwk_ND$colors)

# Plot the dendrogram and the module colors underneath

plotDendroAndColors(netwk_ND$dendrograms[[1]],
                    mergedColors[netwk_ND$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05 )

table(netwk_ND$colors)

# ND eigengenes 
MEs_ND <- moduleEigengenes(datExpr0, colors = 
                             labels2colors(netwk_ND$colors))$eigengenes 
MEs_ND <- orderMEs(MEs_ND) 
colnames(MEs_ND) = names(MEs_ND) %>% gsub("ME","", .)
cor_ND <- cor(MEs_ND)

heatmap.2(cor_ND, 
          main = "ND Module Eigengene Correlation", 
          trace = "none", 
          col = colorRampPalette(c("blue", "white", "red"))(50), # Color gradient 
          key = TRUE,  # Adds a color key (legend) 
          key.title = "Correlation", 
          key.xlab = "Value", 
          density.info = "none",  # Disable histogram in the legend 
          denscol = NA)           # Remove histogram color) 


#Tumor dataset 
TOM_ND <- TOMsimilarityFromExpr(datExpr0, power = 13)

row.names(TOM_ND) <- colnames(datExpr0)
colnames(TOM_ND)  <- colnames(datExpr0)





#---------------------------------------------------------
# 3. Extract phenotype (sample metadata)
#---------------------------------------------------------
pheno <- pData(gse)

# Extract the diabetes status column
pheno$group1 <- pheno[,"diabetes status (nd (non-diabetic), t2d (type 2 diabetic), igt (impaired glucose tolerance, t3cd (type 3 diabetic)):ch1"]

# keep only ND
pheno <- pheno[!pheno$group1 %in% c("ND","IGT","T3cD"), ]

# convert to factor
pheno$group1 <- as.factor(pheno$group1)

# extract selected sample IDs
selected_samples1 <- rownames(pheno)

# Extract expression matrix
#---------------------------------------------------------
expr1 <- exprs(gse)

# keep only selected samples
expr_subset1 <- expr1[, selected_samples1]

dim(expr_subset1)



# 5. Download platform annotation (GPL570)
gpl <- getGEO("GPL570", AnnotGPL = TRUE)
gpl_table <- Table(gpl)


#---------------------------------------------------------
# 6. Prepare expression dataframe for merging
#---------------------------------------------------------
expr_df1 <- as.data.frame(expr_subset1)
expr_df1$ProbeID <- rownames(expr_df1)

# select annotation columns
annotation <- gpl_table[, c("ID", "Gene symbol", "Gene title")]
colnames(annotation) <- c("ProbeID", "GeneSymbol", "GeneTitle")

#---------------------------------------------------------
# 7. Merge probe annotation with expression data
#---------------------------------------------------------
expr_annotated1 <- merge(annotation, expr_df1, by = "ProbeID")


# remove probe column
expr1 <- expr_annotated1[, -1]

# rename gene column
colnames(expr1)[1] <- "Gene"

# remove GeneTitle column
expr1 <- expr1[, -2]


#---------------------------------------------------------
# 8. Collapse duplicate probes by gene (mean expression)
#---------------------------------------------------------
expr_gene1 <- expr1 %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean, na.rm = TRUE))


# 9. Filter out any rows where the Gene name is missing or empty
expr_gene1 <- as.data.frame(expr_gene1)
expr_gene1 <- expr_gene1[!is.na(expr_gene1$Gene) & expr_gene1$Gene != "", ]

# 10. Set the Gene Symbols as the row names
rownames(expr_gene1) <- expr_gene1$Gene
expr_gene1 <- expr_gene1[, -1]
datExpr1 <- t(expr_gene1)

# 2. Check for genes and samples with too many missing values
# WGCNA has a specific function for this cleaning step
library(WGCNA)
gsg <- goodSamplesGenes(datExpr1, verbose = 3)
summary(gsg)
gsg$allOK

# Choose a set of soft-thresholding powers

powers<- c(1:20)
spt1<- pickSoftThreshold(datExpr1, powerVector = powers,
                         verbose = 5,
                         networkType = "signed")
spt1

#plot this data frame to better visualize what β value you should choose.

pdf("Scale_Independence_&_Mean_Connectivity_Plot_T2D.pdf", width = 8, height = 6)
#Plotting the results 

par(mar = c(4.5, 4.5, 2.5, 1.5)) 
#Index the scale-free topology adjustment as a function of the power soft thresholding. 

plot(spt1$fitIndices[,1],-sign(spt1$fitIndices[,3])*spt1$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", 
     main = paste("Scale independence for T2D")) 
text(spt1$fitIndices[,1], -sign(spt1$fitIndices[,3])*spt1$fitIndices[,2], 
     labels=powers,cex=cex1,col="red") 

#This line corresponds to use a cut-off R2 of h 
abline(h=0.8,col="red") 

#Connectivity means as a function of soft power thresholding 
plot(spt1$fitIndices[,1], spt1$fitIndices[,5], 
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", 
     main = paste("Mean connectivity for T2D")) 
text(spt1$fitIndices[,1], spt1$fitIndices[,5], labels=powers, 
     cex=cex1,col="red") 

#This line corresponds to use a cut-off R2 of h 
abline(h=0.8,col="red")

dev.off()

#Calling the Adjacency Function
softPower <- 12
adjacency <- adjacency(datExpr1, power = softPower, type = "signed")

netwk_T2D <- blockwiseModules(
  datExpr1, 
  power = 13,                # The optimal power we found!
  corType = "bicor",         # If using biweight midcorrelation
  maxBlockSize = 10000,
  # Force WGCNA to use its own internal bicor options:
  networkType = "signed",    
  TOMType = "signed",
  verbose = 3
)
mergedColors = labels2colors(netwk_T2D$colors)

# Plot the dendrogram and the module colors underneath

plotDendroAndColors(
  
  netwk_T2D$dendrograms[[1]],
  
  mergedColors[netwk_T2D$blockGenes[[1]]],
  
  "Module colors",
  
  dendroLabels = FALSE,
  
  hang = 0.03,
  
  addGuide = TRUE,
  
  guideHang = 0.05 )

table(netwk_T2D$colors)


# T2D eigengenes 
MEs_T2D <- moduleEigengenes(datExpr1, colors = 
                              labels2colors(netwk_T2D$colors))$eigengenes 
MEs_T2D <- orderMEs(MEs_T2D) 
colnames(MEs_T2D) = names(MEs_T2D) %>% gsub("ME","", .)
cor_T2D <- cor(MEs_T2D)

heatmap.2(cor_T2D, 
          main = "T2D Module Eigengene Correlation", 
          trace = "none", 
          col = colorRampPalette(c("blue", "white", "red"))(50), # Color gradient 
          key = TRUE,  # Adds a color key (legend) 
          key.title = "Correlation", 
          key.xlab = "Value", 
          density.info = "none",  # Disable histogram in the legend 
          denscol = NA)           # Remove histogram color) 


#Tumor dataset 
TOM_T2D <- TOMsimilarityFromExpr(datExpr1, power = 13) 
row.names(TOM_T2D) <- colnames(datExpr1)
colnames(TOM_T2D)  <- colnames(datExpr1)




multiData <- list( 
  T2D = list(data = datExpr1), # Transpose to make samples columns 
  ND = list(data = datExpr0) # Transpose to make samples columns 
) 
# Extract module assignments and colors 

T2D_modules <- netwk_T2D$colors 
T2D_module_colors <- labels2colors(T2D_modules) 
names(T2D_module_colors) <- names(netwk_T2D$colors) 
ND_module_colors <- labels2colors(netwk_ND$colors) 
names(ND_module_colors) <- names(netwk_ND$colors) 

# Module colors for tumor and normal datasets 

multiColor <- list( 
  T2D = T2D_module_colors,   # Named vector of module colors for tumor 
  ND = ND_module_colors  # Named vector of module colors for normal 
) 

#check names to ensure genes name in multiColor match with gene names in the original data 
all(names(multiColor$T2D) %in% rownames(expr_gene1))  # Should return TRUE 

all(names(multiColor$ND) %in% rownames(expr_gene))  # Should return TRUE 


preservation_results <- modulePreservation( 
  multiData = multiData,          
  # List of datasets 
  multiColor = multiColor,        
  referenceNetworks = 1,          
  nPermutations = 300,            
  randomSeed = 12345,             
  verbose = 3                     
) 

#Preservation statistics for modules 
preservation_stats<- 
  preservation_results$preservation$Z$ref.T2D$inColumnsAlsoPresentIn.ND 


#Remove the “gold” module. 
preservation_stats <- preservation_stats[rownames(preservation_stats) != "gold", ]

mod_colors <- rownames(preservation_stats) # Module colors 
Z_summary <- preservation_stats$Zsummary.pres 
preservation_data <- data.frame(Module = mod_colors, Z_summary = Z_summary) 


# Plot Z-summary statistics with point plot 
dev.off()
graphics.off()
ggplot(preservation_data, aes(x = Module, y = Z_summary, color = Module)) + 
  geom_point(size = 8) +  # Larger points 
  scale_color_manual(values = mod_colors) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 1) +  # Threshold line 
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue", size = 1) +  # Strong preservation line 
  theme_minimal() + 
  theme( 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),  # Larger x-axis tick text 
    axis.text.y = element_text(size = 14),  # Larger y-axis tick text 
    axis.title.x = element_text(size = 14, face = "bold"),  # Larger x-axis label 
    axis.title.y = element_text(size = 14, face = "bold"),  # Larger y-axis label 
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  # Title formatting 
  ) + 
  labs( 
    title = "Module Preservation Statistics", 
    x = "Module Colors", 
    y = "Preservation Z-summary" 
  ) 

# Identify highly preserved modules 
highly_preserved <- rownames(preservation_stats)[preservation_stats$Zsummary.pres > 
                                                   10] 
print(highly_preserved)  # List of module colors 

# Identify moderately preserved modules 
moderately_preserved<- 
  rownames(preservation_stats)[preservation_stats$Zsummary.pres >= 3 & preservation_stats$Zsummary.pres <= 10] 
print(moderately_preserved) 

# Identify low-preserved modules 
low_preserved <- rownames(preservation_stats)[preservation_stats$Zsummary.pres < 3] 
print(low_preserved) 


# Run once if extrafont isn't installed
install.packages("extrafont")
library(extrafont)

# Load system fonts (only need to run font_import() once per R installation)
# font_import() 
loadfonts(device = "win") # Use device = "pdf" or device = "postscript" if exporting to PDF


library(ggplot2)
library(dplyr)

# 1. Prepare data frame from preservation statistics
plot_data_all <- data.frame(
  Module = rownames(preservation_stats),
  Z_summary = preservation_stats$Zsummary.pres,
  stringsAsFactors = FALSE
) %>%
  filter(!Module %in% c("gold", "grey")) %>%
  arrange(desc(Z_summary)) %>%
  mutate(Module = factor(Module, levels = Module))

# 2. Build the Main Module Preservation Bar Plot
barplot_Zsummary <- ggplot(plot_data_all, aes(x = Module, y = Z_summary, fill = Module)) +
  geom_col(color = "black", linewidth = 0.3) +
  scale_fill_manual(values = setNames(as.character(plot_data_all$Module), plot_data_all$Module)) +
  
  # Ensure top and bottom elements have room to breathe
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  
  # Preservation threshold lines
  geom_hline(yintercept = 2, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "darkblue", linewidth = 0.8) +
  
  # Threshold Annotations
  annotate("text", x = 0.6, y = 2.8, label = "Low Preservation (Z = 2)", 
           color = "darkred", hjust = 0, fontface = "bold", family = "Times New Roman", size = 3) +
  annotate("text", x = 0.6, y = 10.8, label = "High Preservation (Z = 10)", 
           color = "darkblue", hjust = 0, fontface = "bold", family = "Times New Roman", size = 3) +
  
  theme_minimal(base_family = "Times New Roman") +
  theme(
    # Generous right and top margins so no text/labels are cut off
    plot.margin = margin(t = 20, r = 40, b = 25, l = 20, unit = "pt"),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7.5, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 12, face = "bold"),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 15, vjust = 2),
    plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic", vjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Bar plot of Module Preservation Z-summary Scores",
    subtitle = "Evaluation of network module stability across conditions",
    x = "Modules",
    y = "Preservation Z-summary Score"
  )

# 3. Export to high-resolution JPEG (or PNG)
ggsave(
  filename = "plots_jpeg/Module_Preservation_Barplot.png", 
  plot = barplot_Zsummary, 
  width = 11, 
  height = 6.5, 
  dpi = 300
)





# 1. Prepare plot_data WITH the Alpha column
plot_data <- data.frame(
  Module = rownames(preservation_stats),
  Z_summary = preservation_stats$Zsummary.pres,
  stringsAsFactors = FALSE
) %>%
  filter(!Module %in% c("gold", "grey")) %>%
  mutate(
    # Create the required Alpha column for transparency
    Alpha = ifelse(Module %in% c("darkturquoise", "brown"), 1.0, 0.35)
  ) %>%
  arrange(desc(Z_summary)) %>%
  mutate(Module = factor(Module, levels = Module))
contrast_plot <- ggplot(plot_data, aes(x = Module, y = Z_summary, fill = Module, alpha = Alpha)) +
  geom_col(color = "black", linewidth = 0.3) +
  scale_fill_manual(values = setNames(as.character(plot_data$Module), plot_data$Module)) +
  scale_alpha_identity() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  
  # Preservation lines
  geom_hline(yintercept = 2, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "darkblue", linewidth = 0.8) +
  
  # Brown module score label
  geom_text(
    data = filter(plot_data, Module == "brown"),
    aes(label = sprintf("%.2f", Z_summary)),
    vjust = -0.6, fontface = "bold", family = "Times New Roman", size = 3.5
  ) +
  
  # Darkturquoise score label moved up with a pointer line
  geom_text(
    data = filter(plot_data, Module == "darkturquoise"),
    aes(label = sprintf("Darkturquoise: %.2f", Z_summary)),
    vjust = -3.5, hjust = 0.8, fontface = "bold", color = "darkturquoise", family = "Times New Roman", size = 3.8
  ) +
  annotate(
    "segment", 
    x = which(levels(plot_data$Module) == "darkturquoise"), 
    xend = which(levels(plot_data$Module) == "darkturquoise"),
    y = 8, yend = 2, 
    arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.7
  ) +
  
  theme_minimal(base_family = "Times New Roman") +
  theme(
    # Generous right and top margins so no text is cut off
    plot.margin = margin(t = 20, r = 40, b = 25, l = 20, unit = "pt"),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7.5, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 12, face = "bold"),
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 15, vjust = 2),
    plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic", vjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Module Preservation Contrast Analysis",
    subtitle = "Comparing Darkturquoise (Low Preservation) vs. Brown (High Preservation)",
    x = "WGCNA Modules",
    y = "Preservation Z-summary Score"
  )

ggsave("plots_jpeg/Module_Preservation_Contrast.png", plot = contrast_plot, width = 11, height = 6.5, dpi = 300)

library(knitr)
library(dplyr)

# List of modules specified in Table 10
target_modules <- c(
  "darkturquoise", "brown", "salmon", 
  "darkgreen", "orange", "lightcyan", "midnightblue"
)

# Extract and format the summary table
# Assumes 'preservation_stats' is the output from modulePreservation
table_10 <- preservation_stats %>%
  # Filter for the target modules present in your results
  filter(rownames(.) %in% target_modules) %>%
  mutate(Module = rownames(.)) %>%
  select(Module, Z_summary = Zsummary.pres) %>%
  # Round Z-summary to 2 decimal places
  mutate(Z_summary = round(Z_summary, 2)) %>%
  # Categorize and Interpret based on standard WGCNA thresholds
  mutate(
    `Preservation Category` = case_when(
      Z_summary >= 10 ~ "High",
      Z_summary >= 2  ~ "Moderate",
      Z_summary < 2   ~ "Low / None"
    ),
    Interpretation = case_when(
      Z_summary >= 10 ~ "Strongly preserved network topology across ND and T2D",
      Z_summary >= 2  ~ "Moderately preserved network structure",
      Z_summary < 2   ~ "Non-preserved / T2D-specific network changes"
    )
  ) %>%
  # Reorder columns to match Table 10 format
  select(Module, Z_summary, `Preservation Category`, Interpretation)

# Print as a clean Markdown table for reports/papers
kable(table_10, caption = "Table 10: Module Preservation Scores (ND vs T2D)")

# Create output directory if it doesn't exist
if (!dir.exists("tables")) {
  dir.create("tables", recursive = TRUE)
}

# Export to CSV
write.csv(table_10, file = "tables/Table10_Module_Preservation_Scores.csv", row.names = FALSE)




