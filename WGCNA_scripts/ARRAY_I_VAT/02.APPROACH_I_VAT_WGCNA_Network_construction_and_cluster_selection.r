#===============================================================================
#
#                              Code chunk 1
#                      LOAD DATA AND MAIN DIRECTORY
#
#===============================================================================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "/media/mireia/Disco Duro HDD/Documentos/R/ARRAY_adiposo_6000";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "FINE_01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

#===============================================================================
#
#                              Code chunk 2
#            EXPLORING THE MOST APPROPRIATE SOFT THRESHOLD
#                     BASED ON THE MAIN CONNECTIVITY
#
#===============================================================================
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Adjusting plot axis
axis <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
i = 1
yaxis <- c(axis,1)
xaxis <- c(sft$fitIndices[,1],0)
# Scale-free topology fit index as a function of the soft-thresholding power
plot(xaxis, yaxis,
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# This line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#===============================================================================
#
#                                Code chunk 3
#                  SET ESTABLISH THE POWER OF THE SOFT THRESHOLD
#                 (based on the graph, where the red line cuts off)
#
#===============================================================================
softPower = 16;
adjacency = adjacency(datExpr, power = softPower);

#===============================================================================
#
#                              Code chunk 4
#                     CALCULATE TOM-dissimilarity
#
#===============================================================================
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

#===============================================================================
#
#                                Code chunk 5
#                           PLOT TOM-dissimilarity
#
#===============================================================================
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity Visceral Adipose Tissue (input = 6000 probes)",
     labels = FALSE, hang = 0.04);

#===============================================================================
#
#                                Code chunk 6
#                   CLUSTER/MODULE COLORS AND GENE DENDROGRAM
#
#===============================================================================
# We like large modules, so we set the minimum module size relatively high
minModuleSize = 60;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)

#===============================================================================
#
#                                   Code chunk 7
#                         PLOT DENDROGRAM WITH COLOURS
#                   REPRESENTING CO-EXPRESSION CLUSTERS/MODULES
#
#===============================================================================
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors in Visceral Adipose Tissue (input = 6000 probes)")

#===============================================================================
#
#                            Code chunk 8
#             MERGING MODULES WITH VERY SIMILAR PROFILES
#
#===============================================================================
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes in Visceral Adipose Tissue (input = 6000 probes)",
     xlab = "", sub = "")

#===============================================================================
#
#                                Code chunk 9
#                          CHOOSING FUSION THRESHOLD
#
#===============================================================================
# Merging very similar co-expression clusters/modules
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules/clusters:
mergedMEs = merge$newMEs;

#===============================================================================
#
#                               Code chunk 10
#                               COMPARATIVE OF
#                    DYNAMIC TREE CUT VS MERGED MODULES
#
#===============================================================================
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    cex.rowText = 1.5,
                    cex.colorLabels = 1, cex.dendroLabels = 0.9,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

#===============================================================================
#
#                                Code chunk 11
#                                 SAVE DATA
#
#===============================================================================
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FINE_02-dataInput.RData")
