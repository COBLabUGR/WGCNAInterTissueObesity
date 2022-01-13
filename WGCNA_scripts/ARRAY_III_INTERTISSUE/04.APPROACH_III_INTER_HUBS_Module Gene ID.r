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
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# LOAD OUTPUT OF CODE 3
geneInfo = read.csv2("/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_TISULAR_CONJUNTO/consensusAnalysis-CombinedNetworkResultsWITHGENES.csv", header=TRUE, sep=",", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"));

#===============================================================================
#
#                               Code chunk 2
#                              MODULE GENE ID
#                   (GENES WHOSE MM > 0.8 WILL BE CONSIDERED HUBS)
#
#===============================================================================
COLOR <- names(table(geneInfo$ModuleColor))
nombres <- colnames(geneInfo)
nombres <- gsub("Z_kME_meta_ME","",nombres)
nombres <- gsub("kME_set1_ME","",nombres)
nombres <- gsub("kME_set2_ME","",nombres)
nombres <- gsub("p_","",nombres)
nombres <- gsub("kME_meta_ME","",nombres)

for (i in 1:length(COLOR)) {
color <- COLOR[i]
tosave <- geneInfo[which(geneInfo$ModuleColor == color),c(1:4,grep(paste(color,sep=""), colnames(geneInfo)))]
colorname <- tosave[1,3]
tosave <- geneInfo[which(geneInfo$ModuleColor == color),c(1:4,which(colorname == nombres))]
a <- paste("Z_kME_meta_ME",colorname)
a <- gsub(" ","",a)
MM <- which(colnames(tosave)==a)
tosave <- tosave[order(tosave[,MM]),]

# SAVE DATA
write.csv2(tosave,paste("/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_TISULAR_CONJUNTO/Module Gene ID/",color,"_",colorname,".csv",sep=""))
}
