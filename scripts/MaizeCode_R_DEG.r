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

ref_genes<-read.delim(args[4], header = FALSE, 
                      col.names = c("Chr","Start","Stop","Name","Value","Strand"))
ref_genes<-mutate(ref_genes, GeneID=str_replace(ref_genes$Name, pattern = ".*ID=(gene:)?([^;]+).*", replacement = "\\2")) %>%
  select(-Name, -Value)
ref_genes$GeneID<-str_remove_all(ref_genes$GeneID, pattern = "_.")

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
		select(Chr,Start,Stop,GeneID,logFC,Strand,logCPM,PValue,FDR,Sample,DEG) %>%
		arrange(DEG,Chr,Start)
	write.table(DEGtable,paste0("combined/DEG/DEG_",analysisname,"_",sample1,"_vs_",sample2,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	allDEG<-c(allDEG,DEGtable$GeneID)
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

### To make R object for later plotting gene expression

norm<-cpm(y, normalized.lib.size=T)
genextable<-data.frame(norm, stringsAsFactors = FALSE)
genextable<-mutate(genextable, GID=row.names(genextable))

plot.Expression <- function(gene) {
  
  dataline<-filter(genextable, GID==gene) %>% melt(id=c("GID")) %>%
    select(GID, Replicate=variable, CountPerMillion=value)
  
  uniqueDEGs<-read.delim(paste0("combined/DEG/unique_DEGs_",analysisname,".txt"), header = TRUE)
 
  dataline<-merge(dataline, targets, by=c("Replicate")) %>%
    merge(uniqueDEGs, by=c("GID","Sample"), all.x = TRUE)
  
  dataline$Sample<-as.factor(dataline$Sample)
  dataline$CountPerMillion<-as.numeric(dataline$CountPerMillion)
  dataline<-group_by(dataline, Sample) %>% 
    mutate(Average = mean(CountPerMillion))
  dataline$Average[is.na(dataline$Average)]<-as.numeric(0)
  dataline$DEG<-as.factor(dataline$DEG)
  
  plot<-ggplot(dataline, aes(Sample,DEG)) + 
    geom_col(position="dodge", aes(y=Average, fill=DEG)) + 
    scale_fill_manual(values = c("0"="grey", "UP"="pink", "DOWN"="lightblue"),
                      labels=c("0"="No", "UP"="Up", "DOWN"="Down")) +
    geom_point(aes(y=CountPerMillion), size=2, shape=3) + 
    labs(title = gene, y="cpm") + 
    theme(axis.title.y=element_text(size=10), axis.title.x=element_blank(),
          plot.title=element_text(size=15), 
          axis.text.x=element_text(size=10, angle = 90),
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"))
  plot
  
}

save(plot.Expression,genextable,targets, file = paste0("combined/DEG/ReadyToPlot_",analysisname,".RData"))

