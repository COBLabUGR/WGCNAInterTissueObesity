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
workingDir = "/media/mireia/HardDiskHDD/Documentos/R/ANALISIS_TISULAR_CONJUNTO";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Also load results of network analysis
lnames = load(file = "Consensus-NetworkConstruction-man.RData");
lnames
exprSize = checkSets(multiExpr);
nSets = exprSize$nSets;


#===============================================================================
#
#                              Code chunk 2
#                   QUANTIFY CLUSTER-TRAIT ASSOCIATIONS
#
#===============================================================================
Traits[[1]]$data <- Traits[[1]]$data[,c(-17)]
Traits[[2]]$data <- Traits[[2]]$data[,c(-17)]

# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in 1:nSets)
{
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
#  moduleTraitCor[[set]] <- orderMEs(moduleTraitCor[[set]])
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);

}

nList = length(moduleTraitCor[[1]][,1])

Esq_corr_pos = list()
Esq_corr_neg = list()
Esq_pValue = list()

set=1
for (a in 1:nList)
{
  Corr_pos <- table(moduleTraitCor[[set]][a,4:16] >= 0.8)
  Esq_corr_pos [[a]] <- c(MEColorNames[a],Corr_pos[2])
  Corr_neg <- table(moduleTraitCor[[set]][a,4:16] <= -0.8)
  Esq_corr_neg [[a]] <- c(MEColorNames[a],Corr_neg[2])
  pValue <- table(moduleTraitPvalue[[set]][a,4:16] < 0.1)
  Esq_pValue [[a]] <- c(MEColorNames[a],pValue[2])

}

set=2
for (a in 1:nList)
{
  Corr_pos <- which(moduleTraitCor[[set]][a,4:16] >= 0.8)
  Esq_corr_pos [[a]] <- c(MEColorNames,Corr_pos)
  Corr_neg <- which(moduleTraitCor[[set]][a,4:16] >= 0.8)
  Esq_corr_neg [[a]] <- c(MEColorNames,Corr_neg)
  pValue <- which(moduleTraitCor[[set]][a,4:16] >= 0.8)
  Esq_pValue [[a]] <- c(MEColorNames,pValue)
}



table(moduleTraitCor[[1]][a,4:16] >= 0.8), moduleTraitCor[[1]][a,4:16] <= -0.8)

#===============================================================================
#
#                              Code chunk 3
#    PLOT THE DATA OBTAINED FROM THE ASSOCIATIONS BETWEEN TRAITS AND CLUSTERS
#
#===============================================================================
# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
# Open a suitably sized window (the user should change the window size if necessary)
sizeGrWindow(10,6)
pdf(file = "Plots/ModuleTraitRelationships-female.pdf", wi = 10, he = 7);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                           signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               font.lab.x = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships in", setLabels[set]))
#dev.off();
# Plot the module-trait relationship table for set number 2
set = 2
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                           signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,76
pdf(file = "Plots/ModuleTraitRelationships-male.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               font.lab.x = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships in", setLabels[set]))
#dev.off();

#===============================================================================
#
#                              Code chunk 4
#
#===============================================================================
# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0;
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative]);
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative]);
# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0;
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive]);
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive]);

#===============================================================================
#
#                              Code chunk 5
#
#===============================================================================
textMatrix =  paste(signif(consensusCor, 2), "\n(",
                           signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "Plots/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               font.lab.x = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module-trait relationships across\n",
                            paste(setLabels, collapse = " and ")))

#===============================================================================
#
#                              Code chunk 6
#
#===============================================================================
annot = names(multiExpr[[1]]$data)
# Match probes in the data set to the probe IDs in the annotation file
probes = names(multiExpr[[1]]$data)
probes2annot = match(probes, annot)

#===============================================================================
#
#                                Code chunk 7
#      RELATIONSHIP OF GENES WITH SIGNIFICANT CHARACTERISTICS AND CLUSTERS:
#              GENETIC SIGNIFICANCE (GS) AND MODULE MEMBERSHIP (MM)
#
#===============================================================================
consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list();
kME = list();
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data);
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}

#===============================================================================
#
#                                Code chunk 8
#
#===============================================================================
GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2);
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2);
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);

#===============================================================================
#
#                                Code chunk 9
#
#===============================================================================
GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP);
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 6*nTraits)
rownames(GSmat) = probes;
colnames(GSmat) = spaste(
    c("GS_set1_", "GS_set2_", "p_GS_set1_", "p_GS_set2_", "Z_GS_meta_", "p_GS_meta_"),
    rep(traitNames, rep(6, nTraits)))
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
    c("kME_set1_", "kME_set2_", "p_kME_set1_", "p_kME_set2_", "Z_kME_meta_", "p_kME_meta_"),
    rep(MEnames, rep(6, nMEs)))

GSmat = rbind(GS.metaZ, GS.metaP);
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 2*nTraits)
rownames(GSmat) = probes;
colnames(GSmat) = spaste(
    c("Z_GS_meta_", "p_GS_meta_"),
    rep(traitNames, rep(2, nTraits)))
# Same code for kME:
kMEmat = rbind(kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 2*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
    c("Z_kME_meta_", "p_kME_meta_"),
    rep(MEnames, rep(2, nMEs)))
#===============================================================================
#
#                                Code chunk 10
#                                  SAVE DATA
#
#===============================================================================
info = data.frame(GeneSymbol = annot,
             ModuleLabel = moduleLabels,
             ModuleColor = labels2colors(moduleLabels),
             GSmat,
             kMEmat);
write.csv(info, file = "consensusAnalysis-CombinedNetworkResults_meta.csv",
          row.names = FALSE, quote = FALSE);
