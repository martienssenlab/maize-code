#!/usr/bin/env Rscript

library(ggplot2)
library(ComplexUpset)
library(dplyr)
library(purrr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

degtable<-read.delim(args[1], header = TRUE)
samplename<-args[2]
alltfs<-args[3]
degtable$Group<-factor(degtable$Group, levels=alltfs)
degtable$DEG<-factor(degtable$DEG, levels=rev(c("ns","DOWN","Pairwise_DOWN","Pairwise_Mix","Pairwise_UP","UP")))

tot<-degtable %>% group_by(Group) %>%
  summarize(Tot=sum(Number)) %>%
  mutate(label=paste0("n=",Tot))

plot2<-ggplot(degtable, aes(Group, Number)) +
  geom_col(aes(fill=DEG), position="fill", color="black") +
  scale_fill_manual(values=c("ns"="#2e2e2e","DOWN"="blue","Pairwise_DOWN"="lightblue",
                             "Pairwise_Mix"="grey","Pairwise_UP"="pink", "UP"="red")) +
  geom_text(y=1.04, aes(x=tot$Group, label=tot$label), size=5, data.frame()) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        strip.background = element_rect(fill = 'white', colour = 'white')) +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels = c("0","25%","50%","75%","100%")) +
  coord_cartesian(ylim=c(0,1.04)) +
  labs(fill="", title = "Are the genes bound by each group of TFs DEG in ears?")

pdf(paste0("combined/plots/DEGs_TFs_",samplename,".pdf"),10,8)
print(plot2)
dev.off()
