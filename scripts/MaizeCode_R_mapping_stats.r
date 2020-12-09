#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

summary_stats<-args[1]
analysisname<-args[2]

plot.mapping.stats<-function(stattable, name) {
	  table<-read.delim(stattable, header = TRUE, sep="\t") %>%  
    mutate(Label=paste(Tissue,Sample,Rep)) %>% 
    separate(Passing_filtering, into=c("filt","temp1"), sep=" ") %>%
    separate(All_mapped_reads, into=c("all","temp2"), sep=" ") %>%
    separate(Uniquely_mapped_reads, into=c("unique","temp3"), sep=" ") %>%
    select(-Reference_genome, -temp1, -temp2, -temp3)
  
  table$Line<-as.factor(table$Line)
  table$Tissue<-as.factor(table$Tissue)
  table$Sample<-as.factor(table$Sample)
  table$Rep<-as.factor(table$Rep)
  table$filt<-as.numeric(table$filt)
  table$all<-as.numeric(table$all)
  table$unique<-as.numeric(table$unique)
  
  table<-table %>% mutate(Filtered=Total_reads-filt, unmapped=Total_reads-all, multi=all-unique) %>%
    select(-filt, -all) %>%
    gather(key="Read_type",value="Count", Filtered, unmapped, multi, unique) %>%
    group_by(Line,Tissue,Sample,Rep,Read_type)
  
  table$Count<-as.numeric(table$Count)
  table$Read_type<-factor(table$Read_type, 
                          levels = rev(c("unique","multi","unmapped","Filtered")),
                          labels = rev(c("Uniquely mapped","Multi-mapping","Unmapped","Filtered")))
  table<-arrange(table, desc(Line), desc(Tissue), desc(Sample), desc(Rep))
  order<-unique(table$Label)
  table$Label<-factor(table$Label, levels=rev(order))
  
  plot1<-ggplot(table, aes(Label, Count, fill=Read_type)) +
    geom_bar(stat="identity", position="fill", colour="black", show.legend = T) +
    labs(title="", x="",y="", fill="Reads") +
    scale_fill_manual(values = brewer.pal(4,"Paired")) +
    scale_y_continuous("Percentage", labels=c(0,25,50,75,100)) +
    facet_grid(~Line, scales = "free_x") +
    theme(axis.text.x=element_blank(), 
          panel.grid=element_blank(),
          panel.grid.major.y = element_line(colour="grey"),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size=12))
  
  ymax<-plyr::round_any(max(table$Total_reads), 1e7, f = ceiling) 
  ylims<-c(0, ymax)
  ybreaks<-seq(0, ymax, length.out=5)
  ylabels<-c(round(ybreaks/1e6, 1))
  smax<-max(table$Total_reads)
  smin<-min(table$Total_reads)
  sbreaks<-c(smin,smax)
  slabels<-c(paste(round(smin/1e6, 1), "M in min"),paste(round(smax/1e6, 1), "M in max"))
  
  plot2<-ggplot(table, aes(Label, Count, fill=Read_type)) +
    geom_bar(stat="identity", position="stack", colour="black", show.legend = F) +
    guides(colour=FALSE, fill=FALSE) +
    labs(title = "", x="",y="", fill="") +
    scale_fill_manual(values = brewer.pal(4,"Paired")) +
    scale_y_continuous("Million reads", breaks=ybreaks, limits=ylims, labels=ylabels,
                       sec.axis = dup_axis(name = "", breaks=sbreaks, labels=slabels)) +
    geom_hline(yintercept=smin, color="red", linetype="dashed") + 
    geom_hline(yintercept=smax, color="red", linetype="dashed") +
    facet_grid(~Line, scales = "free_x") +
    theme(axis.text.x=element_text(color="black", size=10, angle=90, vjust=0.5, hjust = 1),
          axis.ticks.x=element_blank(), 
          panel.grid.major.x=element_blank(),
          axis.ticks.y.left=element_line(colour="grey"),
          panel.grid.major.y = element_line(colour="grey"), 
          axis.text.y.right=element_text(colour="red"),
          axis.ticks.y.right=element_line(colour="red"),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  title<-ggdraw() +
    draw_label(paste("Mapping statistics for",name,"samples"), 
               size = 12, fontface="bold")
  
  plot_grid(title, plot1, plot2, ncol=1, align="v", axis="lr", rel_heights = c(0.2,1,2.5))
}  

pdf(paste0("combined/plots/mapping_stats_",analysisname,".pdf"), height=10, width=12)
plot.mapping.stats(summary_stats, analysisname)
dev.off()

