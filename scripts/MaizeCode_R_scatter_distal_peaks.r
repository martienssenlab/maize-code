#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

analysisname<-args[1]
tissue<-args[2]
line<-args[3]
included_samples<-args[4]

table1<-read.delim(args[5], header = FALSE, col.names = c("Tissue","PeakID","PeakQuality","GID","RPKM","Group"))
table1$Group<-factor(table1$Group, levels=c("Top20%","20-40%","40-60%","60-80%","Bottom20%"))
table2<-read.delim(args[6], header = TRUE)
tableTOT<-merge(table1,table2, by="PeakID") %>%
    select(-Chr,Start,Stop)
    
if ( grepl( "RNAseq", included_samples, fixed = TRUE) ) {
  tableTOT<-tableTOT %>%
    mutate(RNAseq=RNAseq_plus+RNAseq_minus) %>%
    rowwise() %>%
    mutate(bothstrandsRNAseq=ifelse(RNAseq_plus>0 && RNAseq_minus>0, "Yes","No")) %>%
    arrange(PeakQuality)
  sum<-group_by(tableTOT,bothstrandsRNAseq) %>% summarize(countRNAseq=n())
  tableTOT<-merge(tableTOT, sum, by="bothstrandsRNAseq", all.x = TRUE) %>%
    rowwise() %>%
    mutate(labelRNAseq=paste0(bothstrandsRNAseq," (",countRNAseq,")"))
  tableTOT$labelRNAseq<-as.factor(tableTOT$labelRNAseq)
  
  plot1<-ggplot(tableTOT, aes(PeakQuality, RNAseq, color=labelRNAseq)) +
    geom_point() +
    scale_color_manual(values=c("grey20","#0394fc")) +
    scale_y_continuous(trans="log10",
      labels=scales::label_number_si(accuracy = 1)) +
    scale_x_continuous(trans="log10",
      labels=scales::label_number_si(accuracy = 1)) +
    labs(title=paste("RNAseq signal at H3K27ac peaks in",line,tissue),
         xaxis="Peak quality (log)",
         yaxis="RNAseq signal (log)",
         color="Is RNAseq signal on both strands?") +
    guides(color=guide_legend(override.aes=list(fill=NA))) +
    theme(title = element_text(size=12),
          axis.title = element_text(size=10),
          panel.grid.major = element_line(color="grey"),
          panel.grid.minor = element_line(color="grey"),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          legend.key=element_blank())
  pdf(paste0("combined/plots/scatter_distal_peaks_RNAseq_vs_PeakQuality_",line,"_",tissue,"_",analysisname,".pdf"),12,8)
  print(plot1)
  dev.off()
}

if ( grepl( "RAMPAGE", included_samples, fixed = TRUE) ) {
  tableTOT<-tableTOT %>%
    mutate(RAMPAGE=RAMPAGE_plus+RAMPAGE_minus) %>%
    rowwise() %>%
    mutate(bothstrandsRAMPAGE=ifelse(RAMPAGE_plus>0 && RAMPAGE_minus>0, "Yes","No")) %>%
    arrange(PeakQuality)
  sum<-group_by(tableTOT,bothstrandsRAMPAGE) %>% summarize(countRAMPAGE=n())
  tableTOT<-merge(tableTOT, sum, by="bothstrandsRAMPAGE", all.x = TRUE) %>%
    rowwise() %>%
    mutate(labelRAMPAGE=paste0(bothstrandsRAMPAGE," (",countRAMPAGE,")"))
  tableTOT$labelRAMPAGE<-as.factor(tableTOT$labelRAMPAGE)
  
  plot2<-ggplot(tableTOT, aes(PeakQuality, RAMPAGE, color=labelRAMPAGE)) +
    geom_point() +
    scale_color_manual(values=c("grey20","#cc00cc")) +
    scale_y_continuous(trans="log10",
      labels=scales::label_number_si(accuracy = 1)) +
    scale_x_continuous(trans="log10",
      labels=scales::label_number_si(accuracy = 1)) +
    labs(title=paste("RAMPAGE signal at H3K27ac peaks in",line,tissue),
         xaxis="Peak quality (log)",
         yaxis="RAMPAGE signal (log)",
         color="Is RAMPAGE signal on both strands?") +
    guides(color=guide_legend(override.aes=list(fill=NA))) +
    theme(title = element_text(size=12),
          axis.title = element_text(size=10),
          panel.grid.major = element_line(color="grey"),
          panel.grid.minor = element_line(color="grey"),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          legend.key=element_blank())
  pdf(paste0("combined/plots/scatter_distal_peaks_RAMPAGE_vs_PeakQuality_",line,"_",tissue,"_",analysisname,".pdf"),12,8)
  print(plot2)
  dev.off()
}

if ( grepl( "shRNA", included_samples, fixed = TRUE) ) {
  tableTOT<-tableTOT %>%
    mutate(shRNA=shRNA_plus+shRNA_minus) %>%
    rowwise() %>%
    mutate(bothstrandsshRNA=ifelse(shRNA_plus>0 && shRNA_minus>0, "Yes","No")) %>%
    arrange(PeakQuality)
  sum<-group_by(tableTOT,bothstrandsshRNA) %>% summarize(countshRNA=n())
  tableTOT<-merge(tableTOT, sum, by="bothstrandsshRNA", all.x = TRUE) %>%
    rowwise() %>%
    mutate(labelshRNA=paste0(bothstrandsshRNA," (",countshRNA,")"))
  tableTOT$labelshRNA<-as.factor(tableTOT$labelshRNA)
  
  plot3<-ggplot(tableTOT, aes(PeakQuality, shRNA, color=labelshRNA)) +
    geom_point() +
    scale_color_manual(values=c("grey20","#cc6600")) +
    scale_y_continuous(trans="log10",
      labels=scales::label_number_si(accuracy = 1)) +
    scale_x_continuous(trans="log10",
      labels=scales::label_number_si(accuracy = 1)) +
    labs(title=paste("shRNA signal at H3K27ac peaks in",line,tissue),
         xaxis="Peak quality (log)",
         yaxis="shRNA signal (log)",
         color="Is shRNA signal on both strands?") +
    guides(color=guide_legend(override.aes=list(fill=NA))) +
    theme(title = element_text(size=12),
          axis.title = element_text(size=10),
          panel.grid.major = element_line(color="grey"),
          panel.grid.minor = element_line(color="grey"),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          legend.key=element_blank())
  pdf(paste0("combined/plots/scatter_distal_peaks_shRNA_vs_PeakQuality_",line,"_",tissue,"_",analysisname,".pdf"),12,8)
  print(plot3)
  dev.off()  
}

if ( grepl( "RNAseq", included_samples, fixed = TRUE) && grepl( "RAMPAGE", included_samples, fixed = TRUE) ) {
  plot4<-ggplot(tableTOT, aes(RNAseq, RAMPAGE, size=PeakQuality, alpha=PeakQuality, color=PeakQuality)) +
  geom_point() +
  scale_color_gradient(low="grey20", high="#cc00cc", guide="legend", trans="log10") +
  scale_alpha_continuous(guide="legend", trans="log10") +
  scale_size_continuous(guide="legend", trans="log10") +
  scale_y_continuous(trans="log10",
                     labels=scales::label_number_si(accuracy = 1)) +
  scale_x_continuous(trans="log10",
                     labels=scales::label_number_si(accuracy = 1)) +
  labs(title=paste("RAMPAGE vs RNAseq signal at H3K27ac peaks in",line,tissue),
       xaxis="RNAseq signal (log)",
       yaxis="RAMPAGE signal (log)",
       color="Peak Quality",
       size="Peak Quality",
       alpha="Peak Quality") +
  theme(title = element_text(size=12),
        axis.title = element_text(size=10),
        panel.grid.major = element_line(color="grey"),
        panel.grid.minor = element_line(color="grey"),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        legend.key=element_blank())
  
  pdf(paste0("combined/plots/scatter_distal_peaks_RAMPAGE_vs_RNAseq_",line,"_",tissue,"_",analysisname,".pdf"),12,8)
  print(plot4)
  dev.off()
}

if ( grepl( "RNAseq", included_samples, fixed = TRUE) && grepl( "shRNA", included_samples, fixed = TRUE) ) {
  plot5<-ggplot(tableTOT, aes(RNAseq, shRNA, size=PeakQuality, alpha=PeakQuality, color=PeakQuality)) +
  geom_point() +
  scale_color_gradient(low="grey20", high="#cc6600", guide="legend", trans="log10") +
  scale_alpha_continuous(guide="legend", trans="log10") +
  scale_size_continuous(guide="legend", trans="log10") +
  scale_y_continuous(trans = "log10",
                     labels=scales::label_number_si(accuracy = 1)) +
  scale_x_continuous(trans = "log10",
                     labels=scales::label_number_si(accuracy = 1)) +
  labs(title=paste("shRNA vs RNAseq signal at H3K27ac peaks in",line,tissue),
       xaxis="RNAseq signal (log)",
       yaxis="shRNA signal (log)",
       color="Peak Quality",
       size="Peak Quality",
       alpha="Peak Quality") +
  theme(title = element_text(size=12),
        axis.title = element_text(size=10),
        panel.grid.major = element_line(color="grey"),
        panel.grid.minor = element_line(color="grey"),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        legend.key=element_blank())
  
  pdf(paste0("combined/plots/scatter_distal_peaks_shRNA_vs_RNAseq_",line,"_",tissue,"_",analysisname,".pdf"),12,8)
  print(plot5)
  dev.off()
}

if ( grepl( "RAMPAGE", included_samples, fixed = TRUE) && grepl( "shRNA", included_samples, fixed = TRUE) ) {
  plot5<-ggplot(tableTOT, aes(RAMPAGE, shRNA, size=PeakQuality, alpha=PeakQuality, color=PeakQuality)) +
  geom_point() +
  scale_color_gradient(low="grey20", high="#0394fc", guide="legend", trans="log10") +
  scale_alpha_continuous(guide="legend", trans="log10") +
  scale_size_continuous(guide="legend", trans="log10") +
  scale_y_continuous(trans="log10",
                     labels=scales::label_number_si(accuracy = 1)) +
  scale_x_continuous(trans="log10",
                     labels=scales::label_number_si(accuracy = 1)) +
  labs(title=paste("shRNA vs RAMPAGE signal at H3K27ac peaks in",line,tissue),
       xaxis="RAMPAGE signal (log)",
       yaxis="shRNA signal (log)",
       color="Peak Quality",
       size="Peak Quality",
       alpha="Peak Quality") +
  theme(title = element_text(size=12),
        axis.title = element_text(size=10),
        panel.grid.major = element_line(color="grey"),
        panel.grid.minor = element_line(color="grey"),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        legend.key=element_blank())
  
  pdf(paste0("combined/plots/scatter_distal_peaks_shRNA_vs_RAMPAGE_",line,"_",tissue,"_",analysisname,".pdf"),12,8)
  print(plot5)
  dev.off()
}
