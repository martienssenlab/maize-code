#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

peak_stats<-args[1]
analysisname<-args[2]

plot.peak.stats<-function(stattable, name) {
  table<-read.delim(stattable, header = TRUE, sep="\t") %>%
    mutate(Label=paste(Line,Tissue,Mark)) %>%
    separate(Common_peaks, into=c("Common","temp1"), sep=" ") %>%
    separate(Common_peaks_IDR_0.05, into=c("IDR","temp2"), sep=" ") %>%
    separate(Selected_peaks, into=c("Selected","temp3"), sep=" ") %>%
    rename(Rep1=Peaks_in_Rep1,Rep2=Peaks_in_Rep2,Merged=Peaks_in_merged, Pseudoreps=Peaks_in_pseudo_reps) %>%
    select(-temp1, -temp2, -temp3)
  table$Line<-as.factor(table$Line)
  table$Tissue<-as.factor(table$Tissue)
  table$Mark<-as.factor(table$Mark)
  table$Common<-as.numeric(table$Common)
  table$IDR<-as.numeric(table$IDR)
  table$Selected<-as.numeric(table$Selected)
  table<-gather(table, key="Peak_type",value="Number",Rep1,Rep2,Common,IDR,Merged,Pseudoreps,Selected)
  table$Peak_type<-factor(table$Peak_type, levels=c("Rep1","Rep2","Common","IDR","Merged","Pseudoreps","Selected"))
  table<-arrange(table, desc(Line), desc(Mark), desc(Tissue))
  
  plot<-ggplot(table, aes(Mark,Number,fill=Peak_type)) +
    geom_bar(stat="identity", position="dodge", color="black", show.legend = T) +
    labs(title=paste("Number of peaks in each ChIPseq sample of",analysisname), 
         x="",y="Number of peaks", fill="Peaks in:") +
    scale_fill_manual(values = brewer.pal(7,"Paired")) +
    facet_grid(~Line+Tissue) +
    theme(axis.text.x=element_text(color="black",size=10, angle=90, vjust=0.5, hjust = 1),
          title = element_text(size=15),
          axis.title.y = element_text(size=12),
          legend.title = element_text(size=12),
          panel.grid=element_blank(),
          panel.grid.major.y = element_line(color="grey"),
          panel.grid.minor.y = element_line(color="grey"),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          strip.text.x = element_text(size=12),
          strip.background = element_blank())

  plot
}


pdf(paste0("combined/plots/peak_stats_",analysisname,".pdf"), height=10, width=12)
plot.peak.stats(peak_stats, analysisname)
dev.off()
