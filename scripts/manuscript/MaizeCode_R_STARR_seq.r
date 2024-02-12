#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(gridExtra)
library(patchwork)

args = commandArgs(trailingOnly=TRUE)

tab<-read.delim(args[1], header = TRUE)
titlelabel<-args[2]
size<-args[3]

tab$Tissue<-as.factor(tab$Tissue)
tab$Group<-factor(tab$Group, levels=c("with_k4_expdouble", "with_k4_expsingle", "with_k4_noexp",
          "without_k4_expdouble", "without_k4_expsingle", "without_k4_noexp","shuffle"))

comps1<-list( c("without_k4_expdouble", "without_k4_expsingle"),
             c("without_k4_expdouble", "without_k4_noexp"),
             c("without_k4_expsingle", "without_k4_noexp") )
			 
comps2<-list( c("without_k4_expdouble", "shuffle"), 
			  c("without_k4_expdouble", "without_k4_expsingle"),
              c("without_k4_expdouble", "without_k4_noexp"),
              c("without_k4_expsingle", "without_k4_noexp") )

colors<-c("with_k4_expdouble" = "#5546EE", "with_k4_expsingle" = "#33A9ED", "with_k4_noexp" = "#33CBED", 
          "without_k4_expdouble" = "#ED4433", "without_k4_expsingle" = "#ED8A33", "without_k4_noexp" = "#EDB533", "shuffle" = "#A9A9A9")

###

plot.max.w<-function(tissue, incl, comps) {
  temp<-filter(tab, str_detect(Group, incl), Tissue==tissue)
  plot0<-ggplot(temp, aes(Group, Max, fill=Group)) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_manual(values = colors) +
      labs(title=paste0(tissue)) +
      theme_bw()
  
  findlimits<-function(x) { boxplot.stats(temp[temp$Group==x,]$Max)$stats[c(1,5)] }
  limits<-c(min(sapply(unique(temp$Group), findlimits)), max(sapply(unique(temp$Group), findlimits)))
  min<-ifelse(limits[1]>0, limits[1]*1.1,limits[1]*0.9)
  max<-ifelse(limits[2]>0, limits[2]*1.1,limits[2]*0.9)
  label_y_positions<-seq(11,7,length.out=length(comps))
  plotfin <- plot0 + coord_cartesian(ylim=c(min,max)) +
    stat_compare_means(comparisons = comps, method = "t.test", paired=FALSE, label="p.signif", label.y = label_y_positions,
                       tip.length = c(0))
  
  print(plotfin)
}

###

plot.mean.w<-function(tissue, incl, comps) {
  temp<-filter(tab, str_detect(Group, incl), Tissue==tissue)
  plot0<-ggplot(temp, aes(Group, Mean, fill=Group)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = colors) +
    labs(title=paste0(tissue)) +
    theme_bw()
  
  findlimits<-function(x) { boxplot.stats(temp[temp$Group==x,]$Mean)$stats[c(1,5)] }
  limits<-c(min(sapply(unique(temp$Group), findlimits)), max(sapply(unique(temp$Group), findlimits)))
  min<-ifelse(limits[1]>0, limits[1]*1.1,limits[1]*0.9)
  max<-ifelse(limits[2]>0, limits[2]*1.1,limits[2]*0.9)
  label_y_positions<-seq(-0.7,-0.8,length.out=length(comps))
  plotfin <- plot0 + coord_cartesian(ylim=c(min,max)) +
    stat_compare_means(comparisons = comps, method = "t.test", paired=FALSE, label="p.signif", label.y = label_y_positions,
                       tip.length = c(0))
  
  print(plotfin)
}

###

plot.median.w<-function(tissue, incl, comps) {
  temp<-filter(tab, str_detect(Group, incl), Tissue==tissue)
  plot0<-ggplot(temp, aes(Group, Median, fill=Group)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = colors) +
    labs(title=paste0(tissue)) +
    theme_bw()
  
  findlimits<-function(x) { boxplot.stats(temp[temp$Group==x,]$Median)$stats[c(1,5)] }
  limits<-c(min(sapply(unique(temp$Group), findlimits)), max(sapply(unique(temp$Group), findlimits)))
  min<-ifelse(limits[1]>0, limits[1]*1.1,limits[1]*0.9)
  max<-ifelse(limits[2]>0, limits[2]*1.1,limits[2]*0.9)
  label_y_positions<-seq(-5,-6,length.out=length(comps))
  plotfin <- plot0 + coord_cartesian(ylim=c(min,max)) +
    stat_compare_means(comparisons = comps, method = "t.test", paired=FALSE, label="p.signif", label.y = label_y_positions,
                       tip.length = c(0))
  
  print(plotfin)
}

## Plot Max + mean + median with and without k4

allwith<-ggarrange(plot.max.w("cn","with",comps1) + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank()),
          plot.max.w("ears","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          plot.max.w("roots","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          plot.max.w("endosperm","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          plot.mean.w("cn","with",comps1) + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank()),
          plot.mean.w("ears","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          plot.mean.w("roots","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          plot.mean.w("endosperm","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          plot.median.w("cn","with",comps1) + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank()),
          plot.median.w("ears","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          plot.median.w("roots","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          plot.median.w("endosperm","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          nrow=3, ncol=4, legend="none")

pdf(paste0("manuscript/plots/enhancers_",titlelabel,"_with.pdf"), paper="a4")
print(allwith)
dev.off()

## Plot Max + mean + meadian without k4 only

allwithout<-ggarrange(plot.max.w("cn","without",comps1) + theme(axis.title.x=element_blank(),
                                                   axis.text.x=element_blank(),
                                                   axis.ticks.x=element_blank()),
                   plot.max.w("ears","without",comps1) + theme(axis.title=element_blank(),
                                                     axis.text.x=element_blank(),
                                                     axis.ticks.x=element_blank()),
                   plot.max.w("roots","without",comps1) + theme(axis.title=element_blank(),
                                                      axis.text.x=element_blank(),
                                                      axis.ticks.x=element_blank()),
                   plot.max.w("endosperm","without",comps1) + theme(axis.title=element_blank(),
                                                          axis.text.x=element_blank(),
                                                          axis.ticks.x=element_blank()),
                   plot.mean.w("cn","without",comps1) + theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_blank(),
                                                    axis.ticks.x=element_blank()),
                   plot.mean.w("ears","without",comps1) + theme(axis.title=element_blank(),
                                                      axis.text.x=element_blank(),
                                                      axis.ticks.x=element_blank()),
                   plot.mean.w("roots","without",comps1) + theme(axis.title=element_blank(),
                                                       axis.text.x=element_blank(),
                                                       axis.ticks.x=element_blank()),
                   plot.mean.w("endosperm","without",comps1) + theme(axis.title=element_blank(),
                                                           axis.text.x=element_blank(),
                                                           axis.ticks.x=element_blank()),
                   plot.median.w("cn","without",comps1) + theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_blank(),
                                                    axis.ticks.x=element_blank()),
                   plot.median.w("ears","without",comps1) + theme(axis.title=element_blank(),
                                                      axis.text.x=element_blank(),
                                                      axis.ticks.x=element_blank()),
                   plot.median.w("roots","without",comps1) + theme(axis.title=element_blank(),
                                                       axis.text.x=element_blank(),
                                                       axis.ticks.x=element_blank()),
                   plot.median.w("endosperm","without",comps1) + theme(axis.title=element_blank(),
                                                           axis.text.x=element_blank(),
                                                           axis.ticks.x=element_blank()),
                   nrow=3, ncol=4, legend="none")
pdf(paste0("manuscript/plots/enhancers_",titlelabel,"_without.pdf"), paper="a4")
print(allwithout)
dev.off()

####

### To plot size (length) of enhancers

if ( length(args) == 3 ) {
	plot.size.w<-function(tissue, incl, comps) {
		temp<-filter(tab, str_detect(Group, incl), Tissue==tissue)
		plot0<-ggplot(temp, aes(Group, Size, fill=Group)) +
		geom_boxplot(outlier.shape = NA) +
		scale_fill_manual(values = colors) +
		labs(title=paste0(tissue)) +
		theme_bw()
  
	findlimits<-function(x) { boxplot.stats(temp[temp$Group==x,]$Size)$stats[c(1,5)] }
	limits<-c(min(sapply(unique(temp$Group), findlimits)), max(sapply(unique(temp$Group), findlimits)))
	min<-ifelse(limits[1]>0, limits[1]*1.1,limits[1]*0.9)
	max<-ifelse(limits[2]>0, limits[2]*1.1,limits[2]*0.9)
	label_y_positions<-seq(3000,2500,length.out=length(comps))
	plotfin <- plot0 + coord_cartesian(ylim=c(min,max)) +
				stat_compare_means(comparisons = comps, method = "t.test", paired=FALSE, label="p.signif", label.y = label_y_positions, tip.length = c(0))
  
	print(plotfin)
	}

## Plot sizes with and without k4

	size<-ggarrange(plot.size.w("cn","with",comps1) + theme(axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
					plot.size.w("ears","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
					plot.size.w("roots","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
					plot.size.w("endosperm","with",comps1) + theme(axis.title=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()),
          nrow=1, ncol=4, legend="none", heights=c(1,1,1))

	pdf(paste0("manuscript/plots/Size_enhancers",titlelabel,".pdf"), height=4, width=10)
	print(size)
	dev.off()

	plot1<-plot.size.w("cn","with|shuffle",comps1) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot2<-plot.size.w("ears","with|shuffle",comps1) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot3<-plot.size.w("endosperm","with|shuffle",comps1) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot4<-plot.size.w("roots","with|shuffle",comps1) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot5<-plot.max.w("cn","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot6<-plot.max.w("ears","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot7<-plot.max.w("endosperm","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot8<-plot.max.w("roots","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot9<-plot.mean.w("cn","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot10<-plot.mean.w("ears","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot11<-plot.mean.w("endosperm","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot12<-plot.mean.w("roots","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot13<-plot.median.w("cn","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot14<-plot.median.w("ears","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot15<-plot.median.w("endosperm","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")
	plot16<-plot.median.w("roots","with|shuffle",comps2) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "none")

	combined_plots <- plot1 + plot2 + plot3 + plot4 +
                 plot5 + plot6 + plot7 + plot8 +
                 plot9 + plot10 + plot11 + plot12 +
                 plot13 + plot14 + plot15 + plot16 +
                 plot_layout(nrow = 4, ncol = 4)
				 
	pdf(paste0("manuscript/plots/enhancers_",titlelabel,"with_size.pdf"), paper="a4")
	print(combined_plots)
	dev.off()
	
}




