#!/usr/bin/env Rscript

library(ggplot2)
library(UpSetR)

args = commandArgs(trailingOnly=TRUE)

inputable<-read.delim(args[1], header = TRUE)
samplename<-args[2]

set1<-colnames(inputable)
set2<-set1[! set1 %in% c("PeakID","Distance")]

plot1<-upset(inputable, point.size = 2, line.size = 1, sets = set2, order.by = "freq", keep.order = TRUE,
      queries = list(list(query=elements, params=list("Distance","Gene_body"), color="red", active=T, 
                          query.name=c("Gene_body"))), query.legend = "top")
						  
pdf(paste0("combined/plots/Upset_",samplename,".pdf"),10,8)
print(plot1)
dev.off()
