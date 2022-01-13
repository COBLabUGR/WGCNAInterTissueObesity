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
workingDir = "/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_TISULAR_CONJUNTO";
setwd(workingDir);
# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv2("/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_TISULAR_CONJUNTO/ARRAY_MUSCULO_6000.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1);
# Read in the male liver data set

maleData = read.csv2("/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_TISULAR_CONJUNTO/ARRAY_ADIPOSO_6000.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1);
# Take a quick look at what is in the data sets (caution, longish output):

femData <- femData[,-c(6:10)]
maleData <- maleData[,-c(1,7:12)]
names(femData) <- c("6","7","8","9","10")
names(maleData) <- c("1","2","3","4","5")
names(femData)
names(maleData)

#===============================================================================
#
#                                 Code chunk 2
#               CREATION OF DATAFRAME WITH GENE EXPRESSION DATA
#                       AND ITS ADJUSTMENT FOR ANALYSIS
#
#===============================================================================
# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Muscle tissue", "Adipose tissue")
shortLabels = c("Muscle", "Adipose")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(femData[])));
names(multiExpr[[1]]$data)
rownames(multiExpr[[1]]$data)
multiExpr[[2]] = list(data = as.data.frame(t(maleData[])));
names(multiExpr[[2]]$data)
rownames(multiExpr[[2]]$data)
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

#===============================================================================
#
#                                Code chunk 3
#                        CHECKING THAT DATA IS CORRECT
#
#===============================================================================
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

#===============================================================================
#
#                                 Code chunk 4
#                                 REMOVE GENES
#
#===============================================================================
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

#===============================================================================
#
#                                 Code chunk 5
#                              SAMPLE CLUSTERING
#
#===============================================================================
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);

# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize

#===============================================================================
#
#                               Code chunk 6
#                        CREATION OF TRAIT VARIABLE
#
#===============================================================================
traitData =  read.csv2("/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_TISULAR_CONJUNTO//bd_fenotipos_obesos_bothtissues.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1);
dim(traitData)
names(traitData)

# Remove columns that hold information we do not need.
allTraits <- traitData

# See how big the traits are and what are the trait and sample names
dim(allTraits)
allTraits <- allTraits [, c(-2,-4,-23,-24)];
names(allTraits)
allTraits$Code

# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$Code);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

#===============================================================================
#
#                               Code chunk 7
#                                SAVE DATA
#
#===============================================================================
save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "Consensus-dataInput.RData");
