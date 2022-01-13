#================================================================#===============================================================================
#
#                              Code chunk 1
#                      LOAD DATA AND MAIN DIRECTORY
#
#===============================================================================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "/media/mireia/HardDiskHDD/Documentos/R/ARRAY_Obesidad/Adiposo";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FINE_01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "FINE_02-dataInput.RData");
lnames

#===============================================================================
#
#                              Code chunk 2
#                   QUANTIFY CLUSTER-TRAIT ASSOCIATIONS
#
#===============================================================================
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#===============================================================================
#
#                              Code chunk 3
#    PLOT THE DATA OBTAINED FROM THE ASSOCIATIONS BETWEEN TRAITS AND CLUSTERS
#
#===============================================================================
sizeGrWindow(50,50)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(7, 14, 3, 4));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               xLabelsAngle = 50,
               zlim = c(-1,1),
               cex.lab.x = 1,
               font.lab.x = 1,
               main = paste("Module-trait relationships in Adipose Obesity (input = 6000 probes)"))

#===============================================================================
#
#                                Code chunk 4
#      RELATIONSHIP OF GENES WITH SIGNIFICANT CHARACTERISTICS AND CLUSTERS:
#              GENETIC SIGNIFICANCE (GS) AND MODULE MEMBERSHIP (MM)
#
#===============================================================================
# Define variable BMI Z-score containing the BMI Z-score column of datTrait
zscore = as.data.frame(datTraits$"BMI.Z.score");
names(zscore) = "zscore"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, zscore, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(zscore), sep="");
names(GSPvalue) = paste("p.GS.", names(zscore), sep="");

#===============================================================================
#
#                                 Code chunk 5
#                    IDENTIFYING GENES WITH HIGH GS AND MM
#
#===============================================================================

module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for BMI adjust by age and weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#===============================================================================
#
#                            Code chunk 6
#         SUMMARY OUTPUT TO ANALYSE ONE OF THE CLUSTERS
#
#===============================================================================
# We have found clusters with high association with our trait of interest, and identified their core players through the Module Membership measure. We now merge this statistical information with the genetic annotation and write a file that summarizes the most important genetic annotation and results. It can be inspected in standard spreadsheet software such as Excel or Open Office Calc. Our expression data are only annotated by probe ID names


# To get the total number of genes present in our total data frame
names(datExpr)
names(datExpr)[moduleColors=="red"]

# It informs us of the genes present in the "brown" cluter and which are present in the merged cluster.
# It will return probe IDs belonging to the brown cluster. To facilitate interpretation of the results, we used a probe annotation file provided by the manufacturer of the expression arrays to connect the probe IDs to universally recognized gene names and Entrez codes.

#===============================================================================
#
#                              Code chunk 8
#          CREATE THE VARIABLE WITH THE NECESSARY INFORMATION THE GENES
#
#===============================================================================
# We now create a data frame that contains the following information for all probes:
# PROBE ID
# GENE SYMBOL
# LOCUS LINK ID (Entrez code),
# CLUSTER NAME
# BMI Z-SCORE AND GENE RELATIONSHIP
# MODULE MEMBERSHIP
# GENETIC SIGNIFICANCE

# Cluester will be ordered by their importance to the BMI Z-score, with the most significant on the left.

probes = names(datExpr)
# Create the starting data frame
geneInfo0 = data.frame(Code = probes,
                      moduleColor = moduleColors,
                      geneTraitSignificance,
                      GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, zscore, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM_", modNames[modOrder[mod]], sep=""),
                       paste("p_MM_", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneInfo <- geneInfo0

rownames(geneInfo) <- gsub("\\/_*","",rownames(geneInfo))

#===============================================================================
#
#                             Code chunk 9
#                               SAVE DATA
#
#===============================================================================
write.csv(geneInfo, file = "geneInfo_VAT_II.csv")
