#!/usr/bin/env Rscript

library(AnnotationForge)
library(rrvgo)
library(topGO)

args = commandArgs(trailingOnly=TRUE)

info<-read.delim(args[1], header=FALSE)
genes<-read.delim(args[2], header=TRUE) %>%
 rowwise() %>%
 mutate(desc=ifelse(Description=="protein_coding",Type,Description),
        typ=ifelse(Description=="protein_coding",Description,Type)) %>%
 select(-Description, -Type) %>%
 rename(Description=desc, Type=typ)

fGOzm<-info[,c(2,5,7)]
colnames(fGOzm)<-c("GID","GO","EVIDENCE")

fSymzm<-select(genes, GID, Type, Description)
fSymzm$ENTREZID <- paste0("ent",fSymzm$GID)

fChrzm<-select(genes, GID, Chr)

makeOrgPackage(gene_info=fSymzm, chromosome=fChrzm, go=fGOzm,
              version="0.1",
              maintainer="user <user@maizecode>",
              author="user <user@maizecode>",
              outputDir = "./combined/GO",
              tax_id = "381124",
              genus = "Zea",
              species = "mays",
              goTable="go")

install.packages("./combined/GO/org.Zmays.eg.db", repos=NULL, type="source")
