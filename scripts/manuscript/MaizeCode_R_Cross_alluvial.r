#!/usr/bin/env Rscript

library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)

table<-read.delim(args[1],header = FALSE , stringsAsFactors = TRUE,
                      col.names = c("GID","Tissue","DEG","Inbred","Inbred_ID")) %>%
  rowwise() %>%
  mutate(strata=paste0(Tissue," ",DEG))
table$strata<-as.factor(table$strata)

table1<-select(table, -strata, -Inbred_ID) %>%
  spread(Inbred, DEG) %>%
  group_by(TIL11, B73, NC350, W22, Tissue) %>%
  summarise(DEGs=n()) %>%
  tibble::rownames_to_column("Groups") %>%
  gather(key="Inbred", value="DEG", TIL11, B73, NC350, W22)

table1$Inbred<-factor(table1$Inbred, levels=c("TIL11","NC350","W22","B73"))
table1$DEG<-factor(table1$DEG, levels=c("UP","homolog_noDEG_UP","Other_DEG","homolog_noDEG_Other_DEG","DOWN","homolog_noDEG_DOWN","no_homolog"))

plot.alluvial.DEG<-function(tissue) {
  
  tab<-filter(table1, Tissue==tissue)
  
  label<-tab %>%
    group_by(Inbred, DEG) %>%
    summarize(s=sum(DEGs)) %>%
    arrange(Inbred, desc(DEG)) %>%
    mutate(Labels=paste0(s))
    
  labels<-label$Labels
  
  plot<-ggplot(tab, aes(x=Inbred, stratum=DEG, alluvium=Groups,
                            fill=DEG, y=DEGs)) +
    geom_flow(aes.flow = "backward") +
    geom_stratum(alpha=0.5, size=0.1) +
#    geom_text(stat = "stratum", size = 2.5, label=labels) +
    scale_fill_manual(values=c("UP"="#D70D26","DOWN"="#121479","Other_DEG"="#56137d",
								"homolog_noDEG_UP"="#d45564",
								"homolog_noDEG_DOWN"="#7A7BBA",
								"homolog_noDEG_Other_DEG"="#994ac7","no_homolog"="grey10"), 
						labels=c("no_homolog"="No homologs","UP"="UP regulated vs all other tissues","homolog_noDEG_UP"="Homolog not DEG",
								"DOWN"="DOWN regulated vs all other tissues","homolog_noDEG_DOWN"="Homolog not DEG",
							   "Other_DEG"="DEG vs at least 1 other tissue","homolog_noDEG_Other_DEG"="Homolog not DEG")) +
    scale_x_discrete(expand = c(0.01, 0.01)) +
    theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.x=element_blank(),
		axis.title.y=element_blank(),
        axis.text.x=element_text(size=12),
		axis.text.y=element_text(size=10),
		plot.title=element_text(hjust = 0.5, size=15, margin=margin(0,0,0,0)),
		legend.title=element_blank(),
		legend.text=element_text(size=12),
		legend.key.size=unit(10, "pt")) +
    ggtitle(toupper(paste0(tissue))) +
	guides(fill = guide_legend(nrow = 3))
  
  plot
}

for (tissue in c("ears","cn","pollen","roots")) {  
  assign(paste0("plot_",tissue), plot.alluvial.DEG(tissue))
}

plot<-ggarrange(plot_pollen, plot_roots, plot_cn, plot_ears, ncol=2, nrow=2, align="hv", common.legend = TRUE, legend="bottom")

pdf("manuscript/plots/Fig5x_Alluvial_TIL11_v2_nolabels.pdf", paper="a4")
print(plot)
dev.off()