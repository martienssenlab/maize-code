#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

summary_stats<-args[1]
analysisname<-args[2]

plot.expression.stats<-function(stattable, name) {
  table<-read.delim(stattable, header = TRUE, sep="\t", check.names = FALSE) %>%
    mutate(Label=paste0(Line,"\n(",Total_annotated_genes," genes)"))
  
  table$Label<-as.factor(table$Label)
  table$Tissue<-as.factor(table$Tissue)
  
  table1<-table %>% select(Label,Tissue, Not_expressed_in_Rep1, `Low_expression_in_Rep1(<1cpm)`, `High_expression_in_Rep1(>1cpm)`) %>%
    rename(No=Not_expressed_in_Rep1, `Low (<1cpm)`=`Low_expression_in_Rep1(<1cpm)`, `High (>1cpm)`=`High_expression_in_Rep1(>1cpm)`) %>%
    gather(key="Type",value="Count", No, `Low (<1cpm)`, `High (>1cpm)`) %>%
    mutate(Sample="Rep1")
  
  table2<-table %>% select(Label,Tissue, Not_expressed_in_Rep2, `Low_expression_in_Rep2(<1cpm)`, `High_expression_in_Rep2(>1cpm)`) %>%
    rename(No=Not_expressed_in_Rep2, `Low (<1cpm)`=`Low_expression_in_Rep2(<1cpm)`, `High (>1cpm)`=`High_expression_in_Rep2(>1cpm)`) %>%
    gather(key="Type",value="Count", No, `Low (<1cpm)`, `High (>1cpm)`) %>%
    mutate(Sample="Rep2")
  
  table3<-table %>% select(Label,Tissue, No_mean, `Low_mean(<1cpm)`, `High_mean(>1cpm)`) %>%
    rename(No=No_mean, `Low (<1cpm)`=`Low_mean(<1cpm)`, `High (>1cpm)`=`High_mean(>1cpm)`) %>%
    gather(key="Type",value="Count", No, `Low (<1cpm)`, `High (>1cpm)`) %>%
    mutate(Sample="Average")
    
  tablef<-rbind(table1,table2,table3)
  
  tablef$Type<-factor(tablef$Type, levels = rev(c("No","Low (<1cpm)","High (>1cpm)")))
  tablef$Sample<-factor(tablef$Sample, levels = c("Rep1","Rep2","Average"))
  tablef<-arrange(tablef, desc(Label), desc(Tissue))

  plot<-ggplot(tablef, aes(Sample, Count, fill=Type)) +
    geom_bar(stat="identity", position="fill", colour="black", show.legend = T) +
    labs(title=paste0("Gene expression levels in RNAseq samples of ", analysisname), 
         x="",y="", fill="Expression levels:") +
    scale_fill_manual(values = brewer.pal(3,"Paired")) +
    scale_y_continuous("Percentage", labels=c(0,25,50,75,100)) +
    facet_grid(~Label+Tissue, scales = "free_x") +
    theme(title = element_text(size=15),
          axis.text.x=element_text(size=10), 
          panel.grid=element_blank(),
          panel.grid.major.y=element_line(colour="grey"),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          strip.background=element_blank(),
          strip.text=element_text(size=12))
  plot
}    

pdf(paste0("combined/plots/gene_expression_stats_",analysisname,".pdf"), height=10, width=12)
plot.expression.stats(summary_stats, analysisname)
dev.off()
