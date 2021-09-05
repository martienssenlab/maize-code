#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

summary_stats<-args[1]
analysisname<-args[2]

plot.shRNA.sizes<-function(stattable) {
	count<-read.delim(stattable, header = TRUE) %>%
		spread(key = Type, value=Count) %>%
		mutate(trim=trimmed-filtered, filt=filtered-mapped) %>%
		select(-trimmed, -filtered) %>%
		rename(trimmed=trim, filtered=filt) %>%
		gather(Type, Count, filtered, trimmed, mapped)
	count$Count<-replace_na(count$Count, 0)
	count$Type<-factor(count$Type, levels=c("trimmed","filtered","mapped"))

	plot <- ggplot(count, aes(Size, Count, fill=Type)) +
			geom_bar(stat="identity", position="stack", color="black", size=0.01) +
			facet_wrap(~Sample, nrow = length(unique(count$Sample)), scales="free_y") +
			scale_fill_manual(labels=c("trimmed"="post-trimming","filtered"="post-filtering","mapped"="mapped"), 
                    values = c("trimmed"="grey","filtered"="blue","mapped"="green")) +
			labs(y="Counts", x="Sizes", fill="") +
			scale_x_continuous(breaks = c(15,18,21,24,35,50,80,120,150)) +
			theme(axis.title.y=element_text(size=15), 
				axis.title.x=element_text(size=15),
				axis.text.x=element_text(size=10),
				panel.grid.major.y = element_line(colour="lightgrey"), 
				panel.grid.minor.y = element_blank(),
				panel.grid.major.x = element_line(colour="lightgrey",size=0.1),
				panel.background = element_rect(fill = "white", colour = "black"),
				strip.background = element_rect(fill = 'white', colour = 'black'),
				legend.key=element_blank())
	plot
}  

pdf(paste0("combined/plots/shRNA_sizes_stats_",analysisname,".pdf"), height=10, width=12)
plot.shRNA.sizes(summary_stats)
dev.off()
