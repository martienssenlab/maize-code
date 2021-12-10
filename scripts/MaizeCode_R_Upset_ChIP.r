#!/usr/bin/env Rscript

library(ggplot2)
library(ComplexUpset)
library(dplyr)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

inputable<-read.delim(args[1], header = TRUE) %>%
  mutate(Distance=Distance+1)
samplename<-args[2]

inputable$Group<-factor(inputable$Group, levels=c("Distal_downstream","Terminator","Gene_body","Promoter","Distal_upstream"))
set1<-colnames(inputable)
sampleCols<-set1[! set1 %in% c("PeakID","Distance","Group")]

H3K27ac<-colnames(inputable[grep(pattern="H3K27ac",x=colnames(inputable))])

if ( length(H3K27ac) > 0) {
  combosK27ac <- map(seq(1:length(H3K27ac)), ~ combn(H3K27ac, ., FUN = c, simplify = FALSE)) %>% 
   unlist(recursive = FALSE)
  queryIntersectionsK27ac <- lmap(combosK27ac, ~ if_else(H3K27ac %in% unlist(.x), .x, list(NA))) %>%
    discard(~ all(is.na(.x)))
  queriesK27ac <- map(queryIntersectionsK27ac, 
                 ~ upset_query(intersect = .x, color = "#EE616E", 
                               fill = "#EE616E", only_components = c('intersections_matrix')))
  queriesK27acset <- map(H3K27ac, ~ upset_query(set = .x, fill = "#EE616E"))
}

H3K4me1<-colnames(inputable[grep(pattern="H3K4me1",x=colnames(inputable))])

if ( length(H3K4me1) > 0) {
  combosK4me1 <- map(seq(1:length(H3K4me1)), ~ combn(H3K4me1, ., FUN = c, simplify = FALSE)) %>% 
    unlist(recursive = FALSE)
  queryIntersectionsK4me1 <- lmap(combosK4me1, ~ if_else(H3K4me1 %in% unlist(.x), .x, list(NA))) %>%
    discard(~ all(is.na(.x)))
  queriesK4me1 <- map(queryIntersectionsK4me1, 
                      ~ upset_query(intersect = .x, color = "#8D9BEE", 
                                   fill = "#8D9BEE", only_components = c('intersections_matrix')))
  queriesK4me1set <- map(H3K4me1, ~ upset_query(set = .x, fill = "#8D9BEE"))
}

H3K4me3<-colnames(inputable[grep(pattern="H3K4me3",x=colnames(inputable))])

if ( length(H3K4me3) > 0) {
  combosK4me3 <- map(seq(1:length(H3K4me3)), ~ combn(H3K4me3, ., FUN = c, simplify = FALSE)) %>% 
    unlist(recursive = FALSE)
  queryIntersectionsK4me3 <- lmap(combosK4me3, ~ if_else(H3K4me3 %in% unlist(.x), .x, list(NA))) %>%
    discard(~ all(is.na(.x)))
  queriesK4me3 <- map(queryIntersectionsK4me3, 
                      ~ upset_query(intersect = .x, color = "#F1C062", 
                                    fill = "#F1C062", only_components = c('intersections_matrix')))
  queriesK4me3set <- map(H3K4me3, ~ upset_query(set = .x, fill = "#F1C062"))
 }

queries<-c(queriesK27ac, queriesK4me1, queriesK4me3)

inputable <- inputable %>% 
  mutate(
    exclusive_mark = case_when(
      rowSums(across(all_of(H3K27ac))) == rowSums(across(all_of(sampleCols))) ~ "H3K27ac",
      rowSums(across(all_of(H3K4me1))) == rowSums(across(all_of(sampleCols))) ~ "H3K4me1",
      rowSums(across(all_of(H3K4me3))) == rowSums(across(all_of(sampleCols))) ~ "H3K4me3",
      TRUE ~ as.character("Mix")
    )
  ) %>% 
  relocate(exclusive_mark, .after = Group)

colmarks<-c("H3K27ac"="#EE616E","H3K4me1"="#8D9BEE","H3K4me3"="#F1C062","Mix"="black")

plot<-upset(inputable, sampleCols, name="Peaks", 
      mode='exclusive_intersection',
      n_intersections=30, 
      sort_sets=FALSE,
      height_ratio = 0.75,
      base_annotations = list(
        'Shared peaks'=intersection_size(
          counts=FALSE, mapping=aes(fill=Group)) +
          scale_fill_manual(values=c("Distal_downstream"="#B8B5B3","Terminator"="#B233FF",
                                     "Gene_body"="#3358FF","Promoter"="#FF33E0","Distal_upstream"="#2e2e2e"),
                            name="Distance category")
      ),
      annotations = list(
        'Distance to closest gene' = (
          ggplot(mapping = aes(x=intersection, y=Distance, fill = exclusive_mark)) +
            geom_violin(scale="width", na.rm=TRUE, color = "black") +
            scale_y_continuous(trans = "log10",
                               labels=scales::label_number_si(accuracy = 1, unit = "bp")) +
            scale_fill_manual(values=colmarks, name="Exclusive marks"))
      ),
      queries = queries,
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

pdf(paste0("combined/plots/Upset_Histone_ChIP_",samplename,".pdf"),10,8)
print(plot)
dev.off()
