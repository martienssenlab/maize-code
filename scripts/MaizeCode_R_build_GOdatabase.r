#!/usr/bin/env Rscript

library(AnnotationForge)
library(rrvgo)
library(topGO)
library(dplyr)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

info<-read.delim(args[1], header=FALSE)
genes<-read.delim(args[2], header=TRUE) %>%
 rowwise() %>%
 mutate(desc=ifelse(Description=="protein_coding",Type,Description),
        typ=ifelse(Description=="protein_coding",Description,Type)) %>%
 select(-Description, -Type) %>%
 rename(Description=desc, Type=typ)

line<-args[3]

if ( grepl("TIL", line, fixed=TRUE) ) { ncbiID=76912
 } else { ncbiID=381124 }

fGOzm<-info[,c(2,5,7)]
colnames(fGOzm)<-c("GID","GO","EVIDENCE")

fSymzm<-select(genes, GID, Type, Description)
fSymzm$ENTREZID <- paste0("ent",fSymzm$GID)

fChrzm<-select(genes, GID, Chr)

makeOrgPackage(gene_info=fSymzm, chromosome=fChrzm, go=fGOzm,
              version="0.1",
              maintainer="user <user@maizecode>",
              author="user <user@maizecode>",
              outputDir = paste0("./combined/GO/",line,"/"),
              tax_id = ncbiID,
              genus = "Zea",
              species = "mays",
              goTable="go")

db<-paste0("./combined/GO/",line)
setwd("db")
install.packages("org.Zmays.eg.db", repos=NULL, type="source")
setwd("../../..")
