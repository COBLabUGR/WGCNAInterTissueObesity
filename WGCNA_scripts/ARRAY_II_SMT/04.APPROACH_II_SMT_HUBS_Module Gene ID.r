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
workingDir = "/media/mireia/Disco Duro HDD/Documentos/R/ARRAY_Obesidad/Musculo";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# LOAD MICROARRAY DATA ALREADY PRE-PROCESSED FROM ITS THE CORRESPONDING REPOSITORY
geneInfo0 = read.csv2("/media/mireia/Disco Duro HDD/Documentos/R/ARRAY_Obesidad/Musculo/ARRAY_MUSCULO_6000.csv", header=TRUE, sep=",", stringsAsFactors=F, dec=".", na.strings=c(""," ","NaN","NA"));
# LOAD OUTPUT OF CODE 3
geneInfo = read.csv2("/media/mireia/Disco Duro HDD/Documentos/R/ARRAY_Obesidad/Musculo/geneInfo_SMT_II.csv", header=TRUE, sep=",", stringsAsFactors=F, dec=".", na.strings=c(""," ","NaN","NA"));

#===============================================================================
#
#                               Code chunk 2
#                              MODULE GENE ID
#                   (GENES WHOSE MM > 0.8 WILL BE CONSIDERED HUBS)
#
#===============================================================================
table(gene0[,2] == geneInfo[,2])
COLOR <- names(table(geneInfo$moduleColor))
nombres <- colnames(geneInfo)
nombres <- gsub("MM_","",nombres)
nombres <- gsub("p_","",nombres)

for (i in 1:length(COLOR)) {
color <- COLOR[i]
tosave <- geneInfo[which(geneInfo$moduleColor == color),c(1:3,which(color == nombres))]
a <- paste("MM_",color)
a <- gsub(" ","",a)
MM <- which(colnames(tosave)==a)
tosave <- tosave[order(-tosave[,MM]),]
colorname <- tosave[1,3]

# SAVE DATA
write.csv2(tosave,paste("/media/mireia/Disco Duro HDD/Documentos/R/ARRAY_Obesidad/Musculo/MODULE_GENE_ID_OBESIDAD_M/",colorname,".csv",sep=""))
}
