#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ComplexUpset)
library(purrr)
library(wesanderson)

args = commandArgs(trailingOnly=TRUE)

analysisname<-args[1]
TELabels<-c(unlist(strsplit(args[2],",")))
AllLabels<-c(TELabels,"Intergenic","Terminator","Gene_body","Promoter")
col<-rev(wes_palette("Cavalcanti1", length(TELabels), type = "continuous"))
col<-c(col, wes_palette("Royal1", 4))
pal<-setNames(col, AllLabels)

### For distribution plot

table<-read.delim(args[3], header = TRUE)
table$Gene<-as.factor(table$Gene)
table$TE<-as.factor(table$TE)
table$Label<-factor(table$Label, levels = AllLabels)
table$Labelcombined<-as.factor(table$Labelcombined)

plot1<-ggplot(table, aes(Tissue, fill=Label)) +
  geom_bar(stat="count", position="stack", show.legend = F) +
  labs(title="", x="",y="Number of shRNA clusters") +
  scale_fill_manual(values = pal) +
  theme(panel.grid=element_blank(),
        panel.grid.major.y = element_line(colour="grey"),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
		axis.text.x=element_blank())

plot2<-ggplot(table, aes(Tissue, fill=Label)) +
  geom_bar(stat="count", position="fill", show.legend = T) +
  labs(title="", x="",y="Percentage of shRNA clusters", fill="Genomic feature") +
  scale_fill_manual(values = pal) +
  theme(panel.grid=element_blank(),
        panel.grid.major.y = element_line(colour="grey"),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(size=15))

title<-ggdraw() +
  draw_label(paste0("Distribution of shRNA clusters in ",analysisname), 
             size = 20, fontface="bold")

pdf(paste0("combined/plots/shRNA_clusters_distribution_",analysisname,".pdf"), height=10, width=12)
plot_grid(title, plot1, plot2, ncol=1, align="v", axis="lr", rel_heights = c(0.2,1,2))
dev.off()

### For Upset plot

inputable<-read.delim(args[4], header = TRUE)
inputable$Label<-factor(inputable$Label, levels = AllLabels)
set1<-colnames(inputable)
sampleCols<-set1[! set1 %in% c("Line","Cluster_ID","Gene","TE","Label","Labelcombined","GID")]

plot<-upset(inputable, sampleCols, name="shRNA clusters", 
      mode='exclusive_intersection',
      n_intersections=10, 
      sort_sets=FALSE,
      height_ratio = 0.75,
      base_annotations = list(
        'Shared shRNA clusters'=intersection_size(
          counts=FALSE, mapping=aes(fill=Label))
	  + scale_fill_manual(values = pal)
      ),
      set_sizes = (upset_set_size() + ylab("Total shRNA clusters") +
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

pdf(paste0("combined/plots/Upset_shRNA_clusters_",analysisname,".pdf"),10,8)
print(plot)
dev.off()

### For violin plot, if expression table is given

if ( length(args) == 6 ) {
	inputable<-read.delim(args[5], header = TRUE)
	inputable$Label<-factor(inputable$Label, levels = AllLabels)
	inputable$GID<-as.factor(inputable$GID)
	inputable$Tissues<-as.factor(inputable$Tissues)

	exptable<-read.delim(args[6], header = TRUE)
	exptable$Tissue<-as.factor(exptable$Tissue)
	exptable$GID<-as.factor(exptable$GID)

	filteredtable<-inputable %>%
	  filter(Label=="Gene_body") %>%
	  select(Tissues,GID) %>%
	  unique() %>%
	  merge(exptable, by="GID") %>%
	  group_by(Tissues)
  
	plot3<-ggplot(filteredtable, aes(Tissue, RPKM+1, fill=Tissue)) +
	  geom_violin(show.legend = T) +
	  scale_y_continuous(trans="log10",
                     labels=scales::label_number(scale_cut = cut_short_scale(), accuracy = 1)) +
	  facet_wrap(~Tissues, scales="free_y") +
	  theme(panel.grid=element_blank(),
	        panel.grid.major.y = element_line(colour="grey"),
	        axis.ticks=element_blank(),
	        panel.background=element_blank(),
	        axis.text.x=element_blank())

	pdf(paste0("combined/plots/RNAseq_expression_in_shRNA_genes_",analysisname,".pdf"),12,10)
	print(plot3)
	dev.off()
}
