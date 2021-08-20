#!/usr/bin/env Rscript

library(ggplot2)
library(ComplexUpset)
library(dplyr)
library(purrr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

h3k27file<-args[1]
samplename<-args[2]

if ( h3k27file == "yes" ) {
	
 inputable<-read.delim(args[3], header = TRUE) %>%
  	mutate(Distance=Distance+1) %>%
  	rowwise() %>%
  	mutate(marked=ifelse(H3K27ac==1, "Yes", "No"))

 inputable$Group<-factor(inputable$Group, levels=c("Distal_downstream","Terminator","Gene_body","Promoter","Distal_upstream"))
 set1<-colnames(inputable)
 sampleCols<-set1[! set1 %in% c("PeakID","Distance","Group","GID","marked")]

 colmarks<-c("Yes"="#EE616E", "No"="#2e2e2e")
	
 combosK27ac <- map(seq(1:length(sampleCols)), ~ combn(sampleCols, ., FUN = c, simplify = FALSE)) %>% 
	unlist(recursive = FALSE)

 combosK27ac<-combosK27ac[str_detect(combosK27ac,pattern="H3K27ac")]

 queriesK27ac <- map(combosK27ac, ~ upset_query(intersect = .x, color = "#EE616E", 
                             fill = "#EE616E", only_components = c('intersections_matrix')))

 plot<-upset(inputable, sampleCols, name="Peaks", 
      mode='exclusive_intersection',
      n_intersections=30,
      height_ratio = 0.75,
      base_annotations = list(
        'Shared peaks'=intersection_size(
          mapping=aes(fill=Group)) +
          scale_fill_manual(values=c("Distal_downstream"="#B8B5B3","Terminator"="#B233FF",
                                     "Gene_body"="#3358FF","Promoter"="#FF33E0","Distal_upstream"="#2e2e2e"),
                            name="Distance category")
      ),
      annotations = list(
        'Distance to closest gene' = (
          ggplot(mapping = aes(x=intersection, y=Distance, fill=marked)) +
            geom_violin(scale="width", na.rm=TRUE, color = "black") +
            scale_y_continuous(trans = "log10",
                               labels=scales::label_number_si(accuracy = 1, unit = "bp")) +
			scale_fill_manual(values=colmarks) + guides(fill = FALSE))),
	  queries = queriesK27ac,
      set_sizes = (upset_set_size() + ylab("Total peaks") +
        theme(axis.text.x = element_text(angle = 45))),
      matrix = (intersection_matrix(geom = geom_point(shape = "circle",size = 3),
          segment = geom_segment(size = 1.5),
          outline_color = list(active = alpha("white", 0),inactive = alpha("white", 0))) +
          scale_color_manual(values = c("TRUE" = "black", "FALSE" = alpha("white", 0)),
            labels = c("TRUE" = "yes", "FALSE" = "no"),
            breaks = c("TRUE", "FALSE"),
            guide = FALSE) +
          theme(axis.ticks = element_blank(),
            panel.grid = element_blank())),
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
} else {

inputable<-read.delim(args[3], header = TRUE) %>%
 	 mutate(Distance=Distance+1)

 inputable$Group<-factor(inputable$Group, levels=c("Distal_downstream","Terminator","Gene_body","Promoter","Distal_upstream"))
 set1<-colnames(inputable)
 sampleCols<-set1[! set1 %in% c("PeakID","Distance","Group","GID")]
	
 plot<-upset(inputable, sampleCols, name="Peaks", 
      mode='exclusive_intersection',
      n_intersections=30,
      height_ratio = 0.75,
      base_annotations = list(
        'Shared peaks'=intersection_size(
          mapping=aes(fill=Group)) +
          scale_fill_manual(values=c("Distal_downstream"="#B8B5B3","Terminator"="#B233FF",
                                     "Gene_body"="#3358FF","Promoter"="#FF33E0","Distal_upstream"="#2e2e2e"),
                            name="Distance category")
      ),
      annotations = list(
        'Distance to closest gene' = (
          ggplot(mapping = aes(x=intersection, y=Distance)) +
            geom_violin(scale="width", na.rm=TRUE, color = "black", fill="#2e2e2e") +
            scale_y_continuous(trans = "log10",
                               labels=scales::label_number_si(accuracy = 1, unit = "bp")))),
      set_sizes = (upset_set_size() + ylab("Total peaks") +
        theme(axis.text.x = element_text(angle = 45))),
      matrix = (intersection_matrix(geom = geom_point(shape = "circle",size = 3),
          segment = geom_segment(size = 1.5),
          outline_color = list(active = alpha("white", 0),inactive = alpha("white", 0))) +
          scale_color_manual(values = c("TRUE" = "black", "FALSE" = alpha("white", 0)),
            labels = c("TRUE" = "yes", "FALSE" = "no"),
            breaks = c("TRUE", "FALSE"),
            guide = FALSE) +
          theme(axis.ticks = element_blank(),
            panel.grid = element_blank())),
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
}

pdf(paste0("combined/plots/Upset_TF_",samplename,".pdf"),15,8)
print(plot)
dev.off()
