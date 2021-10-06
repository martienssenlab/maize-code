#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ComplexUpset)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

analysisname<-args[1]

### For distribution plot
table<-read.delim(args[2], header = TRUE)
table$Gene<-as.factor(table$Gene)
table$TE<-as.factor(table$TE)
table$Label<-factor(table$Label, 
                     levels = c("helitron","LINE_element","LTR_retrotransposon",
                                "SINE_element","solo_LTR","terminal_inverted_repeat_element","Intergenic","Terminator","Gene_body",
                                "Promoter"))
table$Labelcombined<-as.factor(table$Labelcombined)

plot1<-ggplot(table, aes(Tissue, fill=Label)) +
  geom_bar(stat="count", position="stack", show.legend = F) +
  labs(title="", x="",y="Number of RAMPAGE peaks") +
  scale_fill_manual(values=c("Intergenic"="#B8B5B3","Terminator"="#B233FF",
                                     "Gene_body"="#3358FF","Promoter"="#FF33E0","helitron"="#0B6D10","LINE_element"="#B9DCBA","LTR_retrotransposon"="#08AF0F",
                                "SINE_element"="#92EB96","solo_LTR"="#11E119","terminal_inverted_repeat_element"="#184F19"),
                            name="Genomic feature") +
  theme(panel.grid=element_blank(),
        panel.grid.major.y = element_line(colour="grey"),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
		axis.text.x=element_blank())

plot2<-ggplot(table, aes(Tissue, fill=Label)) +
  geom_bar(stat="count", position="fill", show.legend = T) +
  labs(title="", x="",y="Percentage of RAMPAGE peaks", fill="Genomic feature") +
  scale_fill_manual(values=c("Intergenic"="#B8B5B3","Terminator"="#B233FF",
                                     "Gene_body"="#3358FF","Promoter"="#FF33E0","helitron"="#0B6D10","LINE_element"="#B9DCBA","LTR_retrotransposon"="#08AF0F",
                                "SINE_element"="#92EB96","solo_LTR"="#11E119","terminal_inverted_repeat_element"="#184F19"),
                            name="Genomic feature") +
  theme(panel.grid=element_blank(),
        panel.grid.major.y = element_line(colour="grey"),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(size=15))

title<-ggdraw() +
  draw_label(paste0("Distribution of RAMPAGE peaks in ",analysisname), 
             size = 20, fontface="bold")

pdf(paste0("combined/plots/TSS_distribution_",analysisname,".pdf"), height=10, width=12)
plot_grid(title, plot1, plot2, ncol=1, align="v", axis="lr", rel_heights = c(0.2,1,2))
dev.off()

### For Upset plot

inputable<-read.delim(args[3], header = TRUE)
inputable$Label<-factor(inputable$Label, levels = c("helitron","LINE_element","LTR_retrotransposon",
                                "SINE_element","solo_LTR","terminal_inverted_repeat_element","Intergenic","Terminator","Gene_body",
                                "Promoter"))
set1<-colnames(inputable)
sampleCols<-set1[! set1 %in% c("Line","Peak_ID","Gene","TE","Label","Labelcombined")]

plot<-upset(inputable, sampleCols, name="RAMPAGE Peaks", 
      mode='exclusive_intersection',
      n_intersections=10, 
      sort_sets=FALSE,
      height_ratio = 0.75,
      base_annotations = list(
        'Shared TSS'=intersection_size(
          counts=FALSE, mapping=aes(fill=Label)) +
          scale_fill_manual(values=c("Intergenic"="#B8B5B3","Terminator"="#B233FF",
                                     "Gene_body"="#3358FF","Promoter"="#FF33E0","helitron"="#0B6D10","LINE_element"="#B9DCBA","LTR_retrotransposon"="#08AF0F",
                                "SINE_element"="#92EB96","solo_LTR"="#11E119","terminal_inverted_repeat_element"="#184F19"),
                            name="Genomic feature")
      ),
      set_sizes = (upset_set_size() + ylab("Total RAMPAGE Peaks") +
        theme(axis.text.x = element_text(angle = 45))),
      matrix = (intersection_matrix(geom = geom_point(shape = "circle",size = 3),
          segment = geom_segment(size = 1.5),
          outline_color = list(active = alpha("white", 0),inactive = alpha("white", 0))) +
          scale_color_manual(values = c("TRUE" = "black", "FALSE" = alpha("white", 0)),
            labels = c("TRUE" = "yes", "FALSE" = "no"),
            breaks = c("TRUE", "FALSE"),
            guide = "none") +
          theme(axis.ticks = element_blank(),
            panel.grid = element_blank())
      ),
      themes = upset_modify_themes(
        list(
          "default" = theme(
            panel.grid.major.x = element_blank(),
            axis.ticks.y = element_line(size = 0.25, color = "#2e2e2e")
          ),
          "intersections_matrix" = theme(
            panel.grid = element_blank(),
            panel.grid.major.y = element_line(color = c("#CFCCCF", "white"), size = 5)
          ),
          "Intersection size" = theme(
            panel.grid = element_blank(),
          ),
          "overall_sizes" = theme(
            panel.grid = element_blank(),
            axis.ticks.x = element_line(size = 0.25, color = "#2e2e2e")
          )
        )
      ),
      stripes = alpha("white", 0)
)

pdf(paste0("combined/plots/Upset_TSS_",analysisname,".pdf"),10,8)
print(plot)
dev.off()
