#!/usr/bin/env Rscript

library(edgeR)
library(AnnotationForge)
library(rrvgo)
library(dplyr)
library(topGO)
library(purrr)

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
sampletable<-read.delim(args[2], header = FALSE)
samplename<-args[3]

keep.exprs<-rowSums(cpm(genecount)>1)>=2
filtered<-genecount[keep.exprs,]
filtered$GID<-row.names(filtered)

info<-read.delim("/grid/martienssen/data_nlsas/jcahn/Genomes/GO/B73_v4_infoGO.tab", header=FALSE)
fGOzm<-info[,c(2,5,7)]
colnames(fGOzm)<-c("GID","GO","EVIDENCE")
geneid2GO<-fGOzm[,c(1,2)]
rn1<-paste(geneid2GO[,1], sep="")
gene2GO<-geneid2GO[,-1]
names(gene2GO)<-rn1

allGenes<-unique(unlist(filtered$GID))
myInterestedGenes<-unique(unlist(sampletable[,4]))
geneList<-factor(as.integer(allGenes %in% myInterestedGenes))
names(geneList)<-allGenes

getGO<-function(ont) {
  GOdata<-new("topGOdata", 
              ontology = ont, 
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = gene2GO)
  resultFisher<-runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  summary<-GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 1000, numChar=1000)
  return(summary %>% mutate(classicFisher = as.numeric(classicFisher)) %>% filter(classicFisher < 0.01))
}

plotGOs<-function(TopGoResults, ont, name) {
  simMatrix<-calculateSimMatrix(TopGoResults$GO.ID,
                                orgdb="org.Zmays.eg.db",
                                ont=ont,
                                method="Rel")
  if ( !is.null(dim(simMatrix)) ) {
    scores<-setNames(-log10(as.numeric(TopGoResults$classicFisher)), TopGoResults$GO.ID)
    reducedTerms<-reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold = 0.7,
                                  orgdb="org.Zmays.eg.db")
    pdf(paste0("combined/plots/topGO_",ont,"_",name,"_scatter.pdf"), width=8, height=8)
    scatterPlot(simMatrix, reducedTerms, size = "score")
    dev.off()
    pdf(paste0("combined/plots/topGO_",ont,"_",name,"_treemap.pdf"), width=8, height=8)
    treemapPlot(reducedTerms, size = "score")
    dev.off()
  }
}

for ( ont in c("BP","MF","CC") ) {
  TopGOresults<-getGO(ont)
  if ( dim(TopGOresults)[1] > 0 ) {
  plotGOs(TopGOresults, ont, samplename)
    }
}
