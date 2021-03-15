#!/usr/bin/env Rscript

library(edgeR)
library(AnnotationForge)
library(rrvgo)
library(dplyr)
library(topGO)

args = commandArgs(trailingOnly=TRUE)

#### If database as to be built. Will need to be prepared for other people to use.
# info<-read.delim("/grid/martienssen/data_nlsas/jcahn/Genomes/GO/B73_v4_infoGO.tab", header=FALSE)
# genes<-read.delim("B73_genes_info.tab", header=TRUE) %>%
#  rowwise() %>%
#  mutate(desc=ifelse(Description=="protein_coding",Type,Description),
#         typ=ifelse(Description=="protein_coding",Description,Type)) %>%
#  select(-Description, -Type) %>%
#  rename(Description=desc, Type=typ)

# fGOzm<-info[,c(2,5,7)]
# colnames(fGOzm)<-c("GID","GO","EVIDENCE")

# fSymzm<-select(genes, GID, Type, Description)
# fSymzm$ENTREZID <- paste0("ent",fSymzm$GID)

# fChrzm<-select(genes, GID, Chr)

# makeOrgPackage(gene_info=fSymzm, chromosome=fChrzm, go=fGOzm,
#               version="0.1",
#               maintainer="user <user@maizecode>",
#               author="user <user@maizecode>",
#               outputDir = ".",
#               tax_id = "381124",
#               genus = "Zea",
#               species = "mays",
#               goTable="go")

# install.packages("./org.Zmays.eg.db", repos=NULL, type="source")

library(org.Zmays.eg.db)

genecount<-read.delim(args[1], header = TRUE, row.names = "gene_ID")
sampletable<-read.delim(args[2], header = c("Chr","Start","End","GID","Score","Strand"))
samplename<-args[3]

keep.exprs<-rowSums(cpm(genecount)>1)>=2
filtered<-genecount[keep.exprs,]
filtered$GID<-row.names(filtered)

info<-read.delim("/grid/martienssen/data_nlsas/jcahn/Genomes/GO/B73_v4_infoGO.tab", header=FALSE)
fGOzm<-info[,c(2,5,7)]
colnames(fGOzm)<-c("GID","GO","EVIDENCE")
geneid2GO<-fGOzm[,c(1,2)]
rn1<-paste(geneid2GO[,1], sep="")
gene2GO = geneid2GO[,-1]
names(gene2GO)<-rn1

allGenes<-unique(unlist(filtered$GID))
myInterestedGenes<-unique(unlist(sampletable$GID))
geneList<-factor(as.integer(allGenes %in% myInterestedGenes))
names(geneList)<-allGenes

getGO<-function(ont) {
  GOdata<-new("topGOdata", 
              ontology = ont, 
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = gene2GO)
  resultFisher<-runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  summary<-GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher",topNodes = 1000, numChar=1000)
  return(summary %>% mutate(classicFisher = as.numeric(classicFisher)) %>% filter(classicFisher < 0.01))
}

plotGOs<-function(TopGoResults, ont, name) {
  simMatrix<-calculateSimMatrix(TopGoResults$GO.ID,
                                orgdb="org.Zmays.eg.db",
                                ont=ont,
                                method="Rel")
  scores<-setNames(-log10(as.numeric(TopGoResults$classicFisher)), TopGoResults$GO.ID)
  reducedTerms<-reduceSimMatrix(simMatrix,
                                scores,
                                threshold = 0.7,
                                orgdb="org.Zmays.eg.db")
  pdf(paste0("combined/plots/topGO_",ont,"_",name,"_scatter.pdf"), width=8, height=8)
  print(scatterPlot(simMatrix, reducedTerms, size = "score"))
  dev.off()
  pdf(paste0("combined/plots/topGO_",ont,"_",name,"_treemap.pdf"), width=8, height=8)
  treemapPlot(reducedTerms, size = "score")
  dev.off()
}

for ( ont in c("BP","MF","CC") ) {
  TopGOresults<-getGO(ont)
  plotGOs(TopGOresults, ont, samplename)
}
  


targets<-read.delim(args[2], header = TRUE)
samples<-as.factor(targets$Sample)
reps<-as.factor(targets$Replicate)
tissues<-unique(samples)

targets$Color<-as.numeric(targets$Color)

colors<-c("black","blue","red","pink","green","purple","lightblue")
color_samples<-colors[targets$Color]

analysisname<-args[3]

ref_genes<-read.delim(args[4], header = FALSE, 
                      col.names = c("Chr","Start","Stop","Name","Value","Strand"))
gene_names<-row.names(genecount)
ref_genes<-mutate(ref_genes, GeneID=str_extract(ref_genes$Name,gene_names)) %>%
  select(-Name, -Value)

# EdgeR analysis
y<-DGEList(counts=filtered, group = samples)
y<-calcNormFactors(y)

#color_samples<-c()
#for (i in 1:length(tissues)) {
#  color_samples<-c(color_samples, rep(colors[i],2))
#}

pdf(paste0("combined/plots/MDS_",analysisname,"_v1.pdf"),10,8)
plotMDS(y, col=color_samples, labels=samples)
dev.off()

pdf(paste0("combined/plots/MDS_",analysisname,"_v2.pdf"),10,8)
plotMDS(y, col=color_samples, labels=reps)
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
