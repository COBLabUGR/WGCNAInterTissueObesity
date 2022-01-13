#===============================================================================
#
#                              Code chunk 1
#                      LOAD DATA AND MAIN DIRECTORY
#
#===============================================================================
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#BiocManager::install("IlluminaHumanMethylationEPICmanifest")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #Specific for the 850k
library(IlluminaHumanMethylationEPICmanifest)
library(org.Hs.eg.db)
library(GO.db)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #Specific for the 850k
library(IlluminaHumanMethylationEPICmanifest)
#BiocManager::install("reshape2")
#BiocManager::install("Gviz")
#BiocManager::install("missMethyl")
#BiocManager::install("GOplot")
#BiocManager::install("DMRcate")
#BiocManager::install("KEGGREST")
library(limma)
library(RColorBrewer)
library(plyr)
library(ggplot2)
library(reshape2)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #Specific for the 850k
library(IlluminaHumanMethylationEPICmanifest) #Specific for the 850k
library(DMRcate)
library(Gviz)
library(missMethyl)
library(GOplot)
library(RColorBrewer)
library(ggplot2)
library(plyr)
library(readr)
library(KEGGREST)
library(stringr)

# The following packages should be changed according to the affymetrix platform under study
require("hgu133plus2.db")
require("pd.hg.u133.plus.2")

#===============================================================================
#
#                              Code chunk 2
#                   CREATION OF CPGS FILE AND ITS ENTREZID
#
#===============================================================================
.getFlatAnnotation <- function(array.type=c("450K","EPIC"),anno=NULL)
  # flatten 450k or EPIC array annotation
  # Jovana Maksimovic
  # 18 September 2018
  # Updated 18 September 2018
  # Modified version of Belida Phipson's .flattenAnn code
{
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } else {
      anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }
  # Get rid of the non-CpG sites
  ann.keep<-anno[grepl("^cg",anno$Name),]

  # Get rid of CpGs that are not annotated
  missing<-ann.keep$UCSC_RefGene_Name==""
  ann.keep<-ann.keep[!missing,]

  # Get individual gene names for each CpG
  geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
  names(geneslist)<-rownames(ann.keep)
  grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  names(grouplist)<-rownames(ann.keep)
  flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
  flat$symbol<-as.character(flat$symbol)
  flat$group <- as.character(flat$group)
  flat$cpg<- substr(rownames(flat),1,10)
  #flat$cpg <- rownames(flat)
  flat$alias <- suppressWarnings(limma:::alias2SymbolTable(flat$symbol))
  eg <- suppressMessages(select(org.Hs.eg.db, keys=keys(org.Hs.eg.db),
                                columns=c("ENTREZID","SYMBOL"),
                                keytype="ENTREZID"))
  colnames(eg) <- c("gene_id","symbol")
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]

  # Keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
  flat.u
  # This randomly samples only 1 gene ID for multimapping CpGs
  #.reduceMultiMap(flat.u)
}
flat.u <- .getFlatAnnotation("EPIC")
names(flat.u)[5] <- "ENTREZID"

#===============================================================================
#
#                              Code chunk 3
#                   CREATION OF ENTREZ ID AND ITS GO
#
#===============================================================================
  keys <- keys(org.Hs.eg.db, keytype = "ENTREZID")
  GeneID.PathID <- suppressMessages(select(org.Hs.eg.db, keys=keys,
                                           columns=c("ENTREZID","GO","ONTOLOGY","PATH"),
                                           keytype="ENTREZID"))
  d <- !duplicated(GeneID.PathID[, c("ENTREZID", "GO")])
  GeneID.PathID <- GeneID.PathID[d, c(1,2,4,5)]
  GOID.TERM <- suppressMessages(select(GO.db, keys=unique(GeneID.PathID$GO),
                                       columns=c("GOID","ONTOLOGY","TERM"),
                                       keytype="GOID"))
  go <- tapply(GeneID.PathID$ENTREZID, GeneID.PathID$GO, list)
  go <- list(idList=go, idTable=GOID.TERM)
names(go$idTable)[1] <- "GO"
go$idTable <- go$idTable[,-2]

# DISPLAY LIST OF ALL GO TERMS (BP,MF AND CC) AND THEIR ASSOCIATED GENES
str(go$idList)

# DISPLAY LIST OF ALL THE GENES IN THE GENOME AND THEIR ASSOCIATED GO
str(GeneID.PathID)

#===============================================================================
#
#                               Code chunk 4
#                      LOAD THE HUBS FILE OF A GIVEN CLUSTER
#
#===============================================================================

#******************************
#*********** NOTE *************
#******************************
# The cluster name should be changed for each of the clusters obtained from the WGCNA for this approach and tissue.
#------------------------------
#------------------------------
#------------------------------

magenta <- read.csv2("/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_GO/MODULE_GENE_ID_META/magenta.csv", header=TRUE, sep=";", stringsAsFactors=F, dec=",", na.strings=c(""," ","NaN","NA"),row.names=1)


#===============================================================================
#
#                               Code chunk 5
#                             LOAD ANNOTATIONS
#
#===============================================================================
MERGED <- merge.data.frame(GeneID.PathID,ae.annots, by= "ENTREZID")
dim(MERGED)
MERGED_ <- merge.data.frame(MERGED,go$idTable, by= "GO",all.x=TRUE)
dim(MERGED_)
MERGED_ <- MERGED_[order(MERGED_$ONTOLOGY),]
MERGED_
#MERGED_ <- subset(MERGED_, ONTOLOGY == "BP")
x <- org.Hs.egPATH
# Get the entrez gene identifiers that are mapped to a KEGG pathway ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
# Get the PATH for the first five genes
xx[1:5]
# Get the first one
xx[[1]]
}

xx[[which(names(xx) %in% "225")]] #donde el 225 se refiere al identificador entrez id
[1] "02010" "04146"
#indx <- sapply(lst, length)
#indx <- lengths(lst)
#res <- as.data.frame(do.call(rbind,lapply(lst, `length<-`,
#                         max(indx))))
#colnames(res) <- names(lst[[which.max(indx)]])
#res
#names(res) <- paste("KEGG_TERM",1:ncol(res),sep="_")
#res$ENTREZID <- names(xx)
dim(MERGED_)

#===============================================================================
#
#                               Code chunk 6
#                          GO ENRICHMENT ANALYSIS
#
#===============================================================================
interestENTREZ <- ae.annots$ENTREZID[-which(is.na(ae.annots$ENTREZID))]
length(interestENTREZ)
table(interestENTREZ %in% flat.u$ENTREZID)
not <- interestENTREZ[-which(interestENTREZ %in% flat.u$ENTREZID)]
ae.annots[which(ae.annots$ENTREZID %in% not),1]
ae.annots[which(ae.annots$ENTREZID %in% not),3]
cpgs <- flat.u[ which(flat.u$ENTREZID %in% interestENTREZ),]
cpgs <- cpgs[!duplicated(cpgs$symbol),]
cpgs
sort(table(cpgs$symbol))
cpgs$cpg
ann850k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)
dim(ann850k)
str(ann850k)
# Get the significant CpG sites at less than 5% FDR
sigCpGs <- cpgs$cpg
sigCpGs
# Total number of significant CpGs at 5% FDR
length(sigCpGs)
# Get all the CpG sites used in the analysis to form the background
all <- ann850k@listData$Name
# Total number of CpG sites tested
length(all)
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=FALSE,array.type="EPIC",collection="GO")
gst_ <- gst[order(gst$P.DE),]
head(gst_)
TO_SAVE <- gst_[which(!is.na(str_extract(gst_$ONTOLOGY, "BP"))),]
TO_SAVE <- TO_SAVE[order(TO_SAVE$P.DE),]
TO_SAVE[c(1:20),]
TO_SAVES <- TO_SAVE[c(1:25),]
TO_SAVES$GENES_DE <- rep(NA,nrow(TO_SAVES))
go$idList
RUTAS <- rownames(TO_SAVES)
for (i in 1:length(RUTAS)) {
if (RUTAS[i] %in% names(go$idList)) {
ENTREZ_ID <- go$idList[[which(names(go$idList) %in% RUTAS[i])]]
TO_SAVES$GENES_DE[i] <- paste(unique(MERGED_[which(MERGED_$ENTREZID %in% ENTREZ_ID),6]),collapse=",")
}else{}
}
TO_SAVES$GENES_DE

#===============================================================================
#
#                               Code chunk 7
#                         KEGG ENRICHMENT ANALYSIS
#
#===============================================================================
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all,plot.bias=FALSE,array.type="EPIC",collection="KEGG")
gst_ <- gst[order(gst$P.DE),]
head(gst_)
TO_SAVE <- gst_
TO_SAVES <- TO_SAVE[c(1:25),]
TO_SAVES$GENES_DE <- rep(NA,nrow(TO_SAVES))
RUTAS <- rownames(TO_SAVES)
RUTAS <- gsub(".*:","",RUTAS)

for (i in 1:length(RUTAS)) {
INFOKEGG <- keggGet(RUTAS[i])
if ("GENE" %in% names(INFOKEGG[[1]])) {
INFOKEGG[[1]]$GENE
ENTREZ_ID <- INFOKEGG[[1]]$GENE[seq(from=1,to=length(INFOKEGG[[1]]$GENE),by=2)]
GENE_SYMBOLS <- INFOKEGG[[1]]$GENE[seq(from=0,to=length(INFOKEGG[[1]]$GENE),by=2)] #GENE SYMBOLS
GENE_SYMBOLS <-gsub(";.*","",GENE_SYMBOLS)
GENE_SYMBOLS
KEGG_PATH <- data.frame(ENTREZ_ID,GENE_SYMBOLS)
KEGG_PATH[which(KEGG_PATH$ENTREZ_ID %in% MERGED_$ENTREZID),2]
TO_SAVES$GENES_DE[i] <- paste(KEGG_PATH[which(KEGG_PATH$ENTREZ_ID %in% MERGED_$ENTREZID),2],collapse=",")
}else{}
}
TO_SAVES$GENES_DE

#===============================================================================
#
#                             Code chunk 8
#                               SAVE DATA
#
#===============================================================================

#******************************
#*********** NOTE *************
#******************************
# The cluster name should be changed for each of the clusters obtained from the WGCNA for this approach and tissue.
#------------------------------
#------------------------------
#------------------------------

write.csv2(MERGED_, file="/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_GO//ENRICHMENT_RESULTS/ae_annot_simplificado_/ae_annot_simplificado_greenyellow.csv", row.names=TRUE)

write.csv2(TO_SAVES, file="/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_GO/ENRICHMENT_RESULTS/enrichmentGO/greenyellow_enrichmentGO.csv", row.names=TRUE)

write.csv2(TO_SAVES, file="/media/mireia/Disco Duro HDD/Documentos/R/ANALISIS_GO/ENRICHMENT_RESULTS/enrichmentKEGG/greenyellow_enrichmentKEGG.csv", row.names=TRUE)
