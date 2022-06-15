#!/usr/bin/env Rscript

library(AnnotationForge)
library(rrvgo)
library(topGO)
library(purrr)
library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(stringr)
library(gplots)

args = commandArgs(trailingOnly=TRUE)

line<-args[1]
db<-paste0("./combined/GO/",line)
setwd(db)
library(org.Zmays.eg.db)
setwd("../../..")

genecount<-read.delim(args[2], header = TRUE, row.names = "gene_ID")
keep.exprs<-rowSums(cpm(genecount)>1)>=2
filtered<-genecount[keep.exprs,]

targets<-read.delim(args[3], header = TRUE)
samples<-as.factor(targets$Sample)
reps<-as.factor(targets$Replicate)
tissues<-unique(samples)

targets$Color<-as.numeric(targets$Color)

colors<-c("black","blue","red","purple","green","lightblue","grey")
color_samples<-colors[targets$Color]

analysisname<-args[4]

ref_genes<-read.delim(args[5], header = FALSE, 
                      col.names = c("Chr","Start","Stop","Name","Value","Strand"))
ref_genes<-mutate(ref_genes, GeneID=str_replace(ref_genes$Name, pattern = ".*ID=(gene:)?([^;]+).*", replacement = "\\2")) %>%
  select(-Name, -Value)

# EdgeR analysis
y<-DGEList(counts=filtered, group = samples)
y<-calcNormFactors(y)

pdf(paste0("combined/plots/MDS_",analysisname,"_v1.pdf"),10,8)
plotMDS(y, col=color_samples, labels=samples)
dev.off()

pdf(paste0("combined/plots/MDS_",analysisname,"_v2.pdf"),10,8)
plotMDS(y, col=color_samples, labels=reps)
dev.off()

pdf(paste0("combined/plots/MDS_",analysisname,"_v3.pdf"),10,8)
plotMDS(y, col=color_samples, labels=samples, dim.plot=c(2,3))
dev.off()

y<-estimateCommonDisp(y, verbose = TRUE)
y<-estimateTagwiseDisp(y)

pdf(paste0("combined/plots/BCV_",analysisname,".pdf"),10,8)
plotBCV(y)
dev.off()

#### Function to create a fold change tables for all genes
create.FC.table<-function(sample1, sample2, table) {
  et<-exactTest(table, pair = c(sample2,sample1))
  out<-topTags(et, n=Inf, adjust.method="BH")
  table<-as.data.frame(out)
  table$GeneID<-rownames(table)
  table<-mutate(table, Sample=paste0(sample1,"_vs_",sample2))
  table
}

#### Function to create a table of DEGs between a pair of samples
create.DEG.table<-function(sample1, sample2, y) {
  et<-exactTest(y, pair = c(sample2,sample1))
  out<-topTags(et, n=Inf, adjust.method="BH")
  table<-filter(out$table, FDR<=0.05, logFC<=-2 | logFC>=2)
  table$GeneID<-rownames(table)
  table<-mutate(table, DEG=as.factor(ifelse(logFC<0, "DOWN", "UP")),Sample=paste0(sample1,"_vs_",sample2))
  table
}

getGO<-function(ont, name, sampletable) {
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
	unnest(genes), by = "GO.ID")
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
  treemapPlot(reducedTerms, size = "score", force.print.labels = TRUE)
  dev.off()
}

filtered$GID<-row.names(filtered)
allGenes<-unique(unlist(filtered$GID))

info<-read.delim(paste0(db,"/",line,"_infoGO.tab"), header=FALSE)
fGOzm<-info[,c(2,5,7)]
colnames(fGOzm)<-c("GID","GO","EVIDENCE")
geneid2GO<-fGOzm[,c(1,2)]
rn1<-paste(geneid2GO[,1], sep="")
gene2GO<-geneid2GO[,-1]
names(gene2GO)<-rn1

allDEG<-c()
for (i in 1:(length(tissues)-1)) {
  sample1<-tissues[i]
  for (j in (i+1):length(tissues)) {
	sample2<-tissues[j]
	FCtable<-create.FC.table(sample1,sample2,y)
	FCtable<-merge(ref_genes,FCtable,by=c("GeneID")) %>%
		select(Chr,Start,Stop,GeneID,logFC,Strand,logCPM,PValue,FDR,Sample) %>%
		arrange(Chr,Start)
	write.table(FCtable,paste0("combined/DEG/FC_",analysisname,"_",sample1,"_vs_",sample2,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	DEGtable<-create.DEG.table(sample1,sample2,y)
	DEGtable<-merge(ref_genes,DEGtable,by=c("GeneID")) %>%
		select(Chr,Start,Stop,GID=GeneID,logFC,Strand,logCPM,PValue,FDR,Sample,DEG) %>%
		arrange(DEG,Chr,Start)
	write.table(DEGtable,paste0("combined/DEG/DEG_",analysisname,"_",sample1,"_vs_",sample2,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	allDEG<-c(allDEG,DEGtable$GID)
	
	updeg<-filter(DEGtable, DEG=="UP")
	myInterestedGenes<-unique(unlist(updeg$GID))
	geneList<-factor(as.integer(allGenes %in% myInterestedGenes))
	names(geneList)<-allGenes
	samplename<-paste0("UP_in_",sample1,"_vs_",sample2)
	for ( ont in c("BP","MF") ) {
		print(paste0("plotting ",ont," for ",samplename))
		getGO(ont, samplename, updeg)
	} 
	downdeg<-filter(DEGtable, DEG=="DOWN")
	myInterestedGenes<-unique(unlist(downdeg$GID))
	geneList<-factor(as.integer(allGenes %in% myInterestedGenes))
	names(geneList)<-allGenes
	samplename<-paste0("DOWN_in_",sample1,"_vs_",sample2)
	for ( ont in c("BP","MF") ) {
		print(paste0("plotting ",ont," for ",samplename))
		getGO(ont, samplename, downdeg)
	}
  }
}

keepDEG<-unique(allDEG)
logcounts<-cpm(y, log=TRUE)
lcpm<-logcounts[keepDEG,]

pdf(paste0("combined/plots/Heatmap_cpm_",analysisname,".pdf"),10,15)
heatmap.2(lcpm,trace="none",ColSideColors = color_samples,
          main=paste0("Differentially expressed genes in all samples from ",analysisname),
          margins=c(12,2),cexCol=2, labRow = "", col="bluered", srtCol=45,
          lwid=c(1,5),lhei=c(0.5,5,0.1), key.title = "", key.xlab = "log(cpm)")
dev.off()

pdf(paste0("combined/plots/Heatmap_zscore_",analysisname,".pdf"),10,15)
heatmap.2(lcpm,trace="none",ColSideColors = color_samples,
          main=paste0("Differentially expressed genes in all samples from ",analysisname),
          margins=c(12,2),cexCol=2, labRow = "", col="bluered", srtCol=45, scale="row",
          lwid=c(1,5),lhei=c(0.5,5,0.1), key.title = "")
dev.off()

