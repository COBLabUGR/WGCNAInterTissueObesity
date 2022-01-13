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
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv("/media/mireia/Disco Duro HDD/Documentos/R/ARRAY_adiposo_6000/ARRAY_ADIPOSO_6000.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1);
# Take a quick look at what is in the data set:
dim(femData);
names(femData);

#===============================================================================
#
#                                 Code chunk 2
#               CREATION OF DATAFRAME WITH GENE EXPRESSION DATA
#                       AND ITS ADJUSTMENT FOR ANALYSIS
#
#===============================================================================
# Remove the first column containing the names of the probes (no for the numerical calculation))
datExpr0 = as.data.frame(t(femData));
names(datExpr0) = rownames(femData);
rownames(datExpr0) = names(femData);
rownames (datExpr0) = c("1", "2","3","4","5","6","7","8","9","10","11")

#===============================================================================
#
#                                Code chunk 3
#                        CHECKING THAT DATA IS CORRECT
#
#===============================================================================
gsg = goodSamplesGenes(datExpr0, verbose = 3);
# It returns:  Flagging genes and samples with too many missing values...
gsg$allOK
# It should return "TRUE"

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

#===============================================================================
#
#                                 Code chunk 4
#                              SAMPLE CLUSTERING
#
#===============================================================================
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(18,11)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0)
plot(sampleTree, main = "Sample clustering to detect outliers in Visceral Adipose Tissue (input = 6000 probes)", sub="", xlab="",cex.lab = 1.5,cex.axis = 1.5, cex.main = 1.5, )

#===============================================================================
#
#                               Code chunk 5
#                           LOAD OF CLINICAL TRAITS
#
#===============================================================================
datExpr = as.data.frame(datExpr0)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#===============================================================================
#
#                               Code chunk 6
#                        CREATION OF TRAIT VARIABLE
#
#===============================================================================
traitData = read.csv2("bd_fenotipos_adiposo.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"));
dim(traitData)
names(traitData)

# Remove columns that hold information we do not need.
allTraits = traitData
#allTraits = allTraits[, c(1, 3:31) ];

dim(allTraits)
names(allTraits)

# Sorting names of allTraits
allTraits <- allTraits[with(allTraits, order(allTraits$Code)), ]

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Code);
datTraits = allTraits[traitRows, -c(1,2,3,26,27,38,45)];
rownames(datTraits) = allTraits[traitRows, 2];
collectGarbage();

#===============================================================================
#
#                               Code chunk 8
#                                 HEATMAP
#
#===============================================================================
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    colorHeightMax = 0.6,
                    cex.colorLabels = 0.7, cex.dendroLabels = 0.9,
                    cex.rowText = 1.2,
                    textPositions = names(datTraits),
                    main = "Sample dendrogram and trait heatmap in Visceral Adipose Tissue (input = 6000 probes)")

#===============================================================================
#
#                               Code chunk 9
#                                SAVE DATA
#
#===============================================================================
save(datExpr, datTraits, file = "FINE_01-dataInput.RData")
