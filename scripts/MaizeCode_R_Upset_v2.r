#!/usr/bin/env Rscript

library(ggplot2)
library(ComplexUpset)

args = commandArgs(trailingOnly=TRUE)

inputable<-read.delim(args[1], header = TRUE)
samplename<-args[2]

set1<-colnames(inputable)
set2<-set1[! set1 %in% c("PeakID","Distance","Group")]

H3K27ac<-colnames(inputable[grep(pattern="H3K27ac",x=colnames(inputable))])

combosK27ac <- map(seq(1:length(H3K27ac)), ~ combn(H3K27ac, ., FUN = c, simplify = FALSE)) %>% 
  unlist(recursive = FALSE)
queryIntersectionsK27ac <- lmap(combosK27ac, ~ if_else(H3K27ac %in% unlist(.x), .x, list(NA))) %>%
  discard(~ all(is.na(.x)))
queriesK27ac <- map(queryIntersectionsK27ac, 
               ~ upset_query(intersect = .x, color = "red", 
                             fill = "red", only_components = c("intersections_matrix","Distance")))

H3K4me1<-colnames(inputable[grep(pattern="H3K4me1",x=colnames(inputable))])

combosK4me1 <- map(seq(1:length(H3K4me1)), ~ combn(H3K4me1, ., FUN = c, simplify = FALSE)) %>% 
  unlist(recursive = FALSE)
queryIntersectionsK4me1 <- lmap(combosK4me1, ~ if_else(H3K4me1 %in% unlist(.x), .x, list(NA))) %>%
  discard(~ all(is.na(.x)))
queriesK4me1 <- map(queryIntersectionsK4me1, 
                    ~ upset_query(intersect = .x, color = "blue", 
                                  fill = "blue", only_components = c("intersections_matrix","Distance")))

H3K4me3<-colnames(inputable[grep(pattern="H3K4me3",x=colnames(inputable))])

combosK4me3 <- map(seq(1:length(H3K4me3)), ~ combn(H3K4me3, ., FUN = c, simplify = FALSE)) %>% 
  unlist(recursive = FALSE)
queryIntersectionsK4me3 <- lmap(combosK4me3, ~ if_else(H3K4me3 %in% unlist(.x), .x, list(NA))) %>%
  discard(~ all(is.na(.x)))
queriesK4me3 <- map(queryIntersectionsK4me3, 
                    ~ upset_query(intersect = .x, color = "purple", 
                                  fill = "purple", 
                                  only_components = c("intersections_matrix","Distance")))

queries<-c(queriesK27ac, queriesK4me1, queriesK4me3)

plot1<-upset(inputable, set2, name="Peaks",
      mode='exclusive_intersection',
      n_intersections=30,
      sort_sets='FALSE',
      base_annotations = list(
        'Intersection size'=intersection_size(
          counts=FALSE, mapping=aes(fill=Group)
        )),
      annotations = list(
        'Distance'=(
          ggplot(mapping=aes(y=Distance)) +
            geom_boxplot(outlier.alpha = 0.1, na.rm=TRUE)
        )),
      queries = queries
      )


pdf(paste0("combined/plots/Upset_",samplename,"_v2.pdf"),10,8)
print(plot1)
dev.off()
