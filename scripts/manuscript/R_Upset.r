#!/usr/bin/env Rscript

library(ggplot2)
library(ComplexUpset)
library(dplyr)
library(purrr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

inputable<-read.delim(args[1], header = TRUE)
outputname<-args[2]

set1<-colnames(inputable)
set1<-set1[-1]

plot1<-upset(inputable, set1,
      mode='exclusive_intersection',
      n_intersections=30,
      height_ratio = 0.75,
	  keep_empty_groups=TRUE
)

pdf(paste0(outputname,".pdf"),10,8)
print(plot1)
dev.off()