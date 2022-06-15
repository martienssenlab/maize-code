#!/usr/bin/env Rscript

library(edgeR)
library(AnnotationForge)
library(rrvgo)
library(dplyr)
library(topGO)
library(purrr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

line<-args[1]
db<-paste0("./combined/GO/",line)
setwd(db)
library(org.Zmays.eg.db)
setwd("../../..")

genecount<-read.delim(args[2], header = TRUE, row.names = "gene_ID")
sampletable<-read.delim(args[3], header = FALSE)
colnames(sampletable)<-c("Chr","Start","Stop","GID","score","strand")
samplename<-args[4]

keep.exprs<-rowSums(cpm(genecount)>1)>=2
filtered<-genecount[keep.exprs,]
filtered$GID<-row.names(filtered)

info<-read.delim(paste0(db,"/",line,"_infoGO.tab"), header=FALSE)
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
  resultFisherSummary <- summary(attributes(resultFisher)$score <= 0.01)
  nSigTerms<-0
  if (length(resultFisherSummary) == 3) {
  	nSigTerms<-as.integer(resultFisherSummary[[3]])
  }
  summary<-GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nSigTerms, numChar=1000)
  tab<-summary %>%
	rename_with(.cols = starts_with("apply"), .fn = ~ return("classicFisher")) %>%
	mutate(classicFisher = classicFisher %>% str_replace(pattern= "< *1e-30", replacement = "1e-30") %>% as.numeric())
  sigTerms<-tab$GO.ID
  genesInTerms<-genesInTerm(GOdata, sigTerms)
  genesInTerms2<-map(genesInTerms, ~ intersect(.x, myInterestedGenes) %>% paste(collapse = ","))
  tab2<-tab %>% 
	left_join(tibble(GO.ID = names(genesInTerms2), genes = genesInTerms2) %>% 
	tidyr::unnest(genes), by = "GO.ID")
  tab3<-tab %>%
	  rename(GO=GO.ID) %>%
	  merge(geneid2GO, by="GO") %>%
	  merge(sampletable, by="GID") %>%
	  select(Chr, Start, Stop, GID, GO, Term) %>%
	  arrange(GO) %>%
	  unique()
  if (nrow(tab2) > 1) {
	  write.table(tab2,paste0(db,"/topGO_",name,"_",ont,"_GOs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  }
  if (nrow(tab3) > 0) {
	  write.table(tab3,paste0(db,"/topGO_",name,"_",ont,"_GIDs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  }	
  
  scores<-setNames(-log10(as.numeric(tab$classicFisher)), tab$GO.ID)
  reducedTerms <- tab2 %>%
          rename("go" = "GO.ID", "term" = "Term") %>%
          mutate(parentTerm = term, score = scores)
  
  if ( nrow(tab) > 1 ) {
  	simMatrix<-calculateSimMatrix(tab$GO.ID,
                                orgdb="org.Zmays.eg.db",
                                ont=ont,
                                method="Rel")
	if ( nrow(simMatrix) > 0 ) {
  	reducedTerms<-reduceSimMatrix(simMatrix,
  	                              scores,
  	                              threshold = 0.7,
  	                              orgdb="org.Zmays.eg.db")
	}
  }
  pdf(paste0("combined/plots/topGO_",name,"_",ont,"_treemap.pdf"), width=8, height=8)
  treemapPlot(reducedTerms, size = "score")
  dev.off()
}

for ( ont in c("BP","MF") ) {
  print(paste0("Getting ",ont," for ",samplename))
  getGO(ont, samplename)
}
