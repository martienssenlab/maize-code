#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

analysisname<-args[1]
table<-read.delim(args[2], header = TRUE)
table$Gene<-as.factor(table$Gene)
table$TE<-as.factor(table$TE)
table$Label<-factor(table$Label, 
                     levels = c("helitron","LINE_element","LTR_retrotransposon",
                                "SINE_element","solo_LTR","terminal_inverted_repeat_element","Intergenic","Terminator","Gene_body",
                                "Promoter"))
table$Labelcombined<-as.factor(table$Labelcombined)

plot1<-ggplot(table, aes(Tissue, fill=Label)) +
  geom_bar(stat="count", position="stack", colour="black", show.legend = F) +
  labs(title="", x="",y="Number of peaks") +
  scale_fill_manual(values = brewer.pal(10,"Paired")) +
  theme(panel.grid=element_blank(),
        panel.grid.major.y = element_line(colour="grey"),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
		axis.text.x=element_blank())
plot1

plot2<-ggplot(table, aes(Tissue, fill=Label)) +
  geom_bar(stat="count", position="fill", colour="black", show.legend = T) +
  labs(title="", x="",y="Percentage", fill="Genomic feature") +
  scale_fill_manual(values = brewer.pal(10,"Paired")) +
  theme(panel.grid=element_blank(),
        panel.grid.major.y = element_line(colour="grey"),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(size=15))
plot2

title<-ggdraw() +
  draw_label(paste0("Distribution of RAMPAGE peaks in ",analysisname), 
             size = 20, fontface="bold")

pdf(paste0("combined/plots/TSS_distribution_",analysisname,".pdf"), height=10, width=12)
plot_grid(title, plot1, plot2, ncol=1, align="v", axis="lr", rel_heights = c(0.2,1,2))
dev.off()
