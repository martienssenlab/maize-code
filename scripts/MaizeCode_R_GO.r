#!/usr/bin/env Rscript

library(edgeR)
library(AnnotationForge)
library(rrvgo)
library(dplyr)
library(topGO)
library(purrr)
library(org.Zmays.eg.db)

args = commandArgs(trailingOnly=TRUE)

genecount<-read.delim(args[1], header = TRUE, row.names = "gene_ID")
sampletable<-read.delim(args[2], header = FALSE)
colnames(sampletable)<-c("Chr","Start","Stop","GID","score","strand")
samplename<-args[3]

keep.exprs<-rowSums(cpm(genecount)>1)>=2
filtered<-genecount[keep.exprs,]
filtered$GID<-row.names(filtered)

info<-read.delim("combined/GO/B73_v4_infoGO.tab", header=FALSE)
fGOzm<-info[,c(2,5,7)]
colnames(fGOzm)<-c("GID","GO","EVIDENCE")
geneid2GO<-fGOzm[,c(1,2)]
rn1<-paste(geneid2GO[,1], sep="")
gene2GO<-geneid2GO[,-1]
names(gene2GO)<-rn1

allGenes<-unique(unlist(filtered$GID))
myInterestedGenes<-unique(unlist(sampletable$GID))
geneList<-factor(as.integer(allGenes %in% myInterestedGenes))
names(geneList)<-allGenes

getGO<-function(ont, name) {
  GOdata<-new("topGOdata", 
              ontology = ont, 
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = gene2GO)
  resultFisher<-runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  nodenb<-min(1000, length(resultFisher@score))
  summary<-GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nodenb, numChar=1000)
  tab<-summary %>%
  	mutate(classicFisher = as.numeric(classicFisher)) %>%
  	filter(classicFisher < 0.01) 
  tab2<-tab %>%
	  rename(GO=GO.ID) %>%
	  merge(geneid2GO, by="GO") %>%
	  merge(sampletable, by="GID") %>%
	  select(Chr, Start, Stop, GID, GO, Term) %>%
	  arrange(GO) %>%
	  unique()
  if (nrow(tab2) > 0) {
	  write.table(tab2,paste0("combined/DEG/topGO_",name,"_",ont,"_GIDs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  }  
  return(tab)
}

plotGOs<-function(TopGoResults, ont, name) {
  simMatrix<-calculateSimMatrix(TopGoResults$GO.ID,
                                orgdb="org.Zmays.eg.db",
                                ont=ont,
                                method="Rel")
  if ( nrow(TopGOresults) > 1 ) {
  	scores<-setNames(-log10(as.numeric(TopGoResults$classicFisher)), TopGoResults$GO.ID)
  	reducedTerms<-reduceSimMatrix(simMatrix,
  	                              scores,
  	                              threshold = 0.7,
  	                              orgdb="org.Zmays.eg.db")
  	# pdf(paste0("combined/plots/topGO_",ont,"_",name,"_scatter.pdf"), width=8, height=8)
  	# scatterPlot(simMatrix, reducedTerms, size = "score")
  	# dev.off()
  	pdf(paste0("combined/plots/topGO_",name,"_",ont,"_treemap.pdf"), width=8, height=8)
  	treemapPlot(reducedTerms, size = "score")
  	dev.off()
  }
}

for ( ont in c("BP","MF") ) {
  TopGOresults<-getGO(ont, samplename)
  if ( nrow(TopGOresults) > 1 ) {
  plotGOs(TopGOresults, ont, samplename)
    }
}
