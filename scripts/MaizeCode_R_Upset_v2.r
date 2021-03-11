#!/usr/bin/env Rscript

library(ggplot2)
library(ComplexUpset)

args = commandArgs(trailingOnly=TRUE)

inputable<-read.delim(args[1], header = TRUE)
samplename<-args[2]

set1<-colnames(inputable)
set2<-set1[! set1 %in% c("PeakID","Distance","Group")]

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
        ))
      )


pdf(paste0("combined/plots/Upset_",samplename,"_v2.pdf"),10,8)
print(plot1)
dev.off()
