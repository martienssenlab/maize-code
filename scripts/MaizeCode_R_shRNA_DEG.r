#!/usr/bin/env Rscript

library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(stringr)
library(gplots)

args = commandArgs(trailingOnly=TRUE)

genecount<-read.delim(args[1], header = TRUE, row.names = "gene_ID")
keep.exprs<-rowSums(cpm(genecount)>1)>=2
filtered<-genecount[keep.exprs,]

targets<-read.delim(args[2], header = TRUE)
samples<-as.factor(targets$Sample)
reps<-as.factor(targets$Replicate)
tissues<-unique(samples)

targets$Color<-as.numeric(targets$Color)

colors<-c("black","blue","red","purple","green","lightblue","grey")
color_samples<-colors[targets$Color]

analysisname<-args[3]

# EdgeR analysis
y<-DGEList(counts=filtered, group = samples)
y<-calcNormFactors(y)

pdf(paste0("plots/MDS_",analysisname,"_v1.pdf"),10,8)
plotMDS(y, col=color_samples, labels=samples)
dev.off()

pdf(paste0("plots/MDS_",analysisname,"_v2.pdf"),10,8)
plotMDS(y, col=color_samples, labels=reps)
dev.off()

pdf(paste0("plots/MDS_",analysisname,"_v3.pdf"),10,8)
plotMDS(y, col=color_samples, labels=samples, dim.plot=c(2,3))
dev.off()

y<-estimateCommonDisp(y, verbose = TRUE)
y<-estimateTagwiseDisp(y)

pdf(paste0("plots/BCV_",analysisname,".pdf"),10,8)
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

allDEG<-c()
for (i in 1:(length(tissues)-1)) {
  sample1<-tissues[i]
  for (j in (i+1):length(tissues)) {
	sample2<-tissues[j]
	FCtable<-create.FC.table(sample1,sample2,y)
	write.table(FCtable,paste0("DEG/FC_",analysisname,"_",sample1,"_vs_",sample2,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	DEGtable<-create.DEG.table(sample1,sample2,y)
	write.table(DEGtable,paste0("DEG/DEG_",analysisname,"_",sample1,"_vs_",sample2,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	allDEG<-c(allDEG,DEGtable$GeneID)
  }
}

keepDEG<-unique(allDEG)
logcounts<-cpm(y, log=TRUE)
lcpm<-logcounts[keepDEG,]

if ( nrow(lcpm) < 65536 ) {

	pdf(paste0("plots/Heatmap_cpm_",analysisname,".pdf"),10,15)
	heatmap.2(lcpm,trace="none",ColSideColors = color_samples,
          main=paste0("Differentially expressed shRNA loci in all samples from ",analysisname),
          margins=c(12,2),cexCol=2, labRow = "", col="bluered", srtCol=45,
          lwid=c(1,5),lhei=c(0.5,5,0.1), key.title = "", key.xlab = "log(cpm)")
	dev.off()

	pdf(paste0("plots/Heatmap_zscore_",analysisname,".pdf"),10,15)
	heatmap.2(lcpm,trace="none",ColSideColors = color_samples,
          main=paste0("Differentially expressed shRNA loci in all samples from ",analysisname),
          margins=c(12,2),cexCol=2, labRow = "", col="bluered", srtCol=45, scale="row",
          lwid=c(1,5),lhei=c(0.5,5,0.1), key.title = "")
	dev.off()
	
}
