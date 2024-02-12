#!/usr/bin/env Rscript

library(Gviz)
library(GenomicFeatures)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)

filenames<-read.delim(args[1], header=TRUE)
if ( file.exists(args[2]) ) {
	genes<-makeTxDbFromGFF(args[2], format="gff3")
} else {
	genes<-c()
}

tes<-import(args[3])
id<-args[4]

htcol<-c()
htcol2<-c()
if ( length(args) == 6 ) {
	htstarttable<-read.delim(args[5], header=FALSE)
	htwidthtable<-read.delim(args[6], header=FALSE)
	colors<-c("#B7E2FD","#fac0c7","#fac0c7","#fac0c7","#fac0c7","#fac0c7")
	colors2<-c("#F6FBFE","#fffafa","#fffafa","#fffafa","#fffafa","#fffafa")
	htstart<-c()
	htwidth<-c()
	for ( i in c(1:nrow(htstarttable)) ) {
		htstart<-c(htstart,htstarttable[i,])
		htwidth<-c(htwidth,htwidthtable[i,])
		htcol<-c(htcol,colors[i])
		htcol2<-c(htcol2,colors2[i])
	}
}

options(ucscChromosomeNames=FALSE)

tracksize<-c(1,1,0.5)
tracklist<-list()
for ( i in c(1:nrow(filenames)) ) {
	name<-filenames$Name[i]
	group<-filenames$Group[i]
	backcolor<-filenames$Backcolor[i]
	trackcolor<-filenames$Trackcolor[i]
	fillcolor1<-filenames$Fillcolor1[i]
	fillcolor2<-filenames$Fillcolor2[i]
	ymin<-filenames$Ymin[i]
	ymax<-filenames$Ymax[i]
	ymintick<-sign(ymin)*((floor(abs(ymin)*100)/100))
	ymaxtick<-sign(ymax)*((floor(abs(ymax)*100)/100))
	tracksize<-c(tracksize,1)
	sample<-paste0(name,"_",group)
	filename<-paste0("data/",sample,"_locus.bw")
	print(paste0("Importing bw for ",sample))
	bw<-import(filename)
	print(paste0("Creating track for ",sample))
	track<-DataTrack(bw, type="polygon", baseline=0, name=sample, background.title = backcolor, col=trackcolor, fill.mountain=c(fillcolor1,fillcolor2), col.baseline="grey50", ylim=c(ymin,ymax), yTicksAt=c(ymintick,ymaxtick), rotation.title=0, cex.title=0.5, lwd=0.01, hjust.title=1)
	tracklist<-append(tracklist, track)
}

axistrack<-GenomeAxisTrack(scale=0.1, labelPos = "above")
genetrack<-GeneRegionTrack(genes, name = "Genes", shape="smallArrow", col="black", fill="grey60", rotation.title=0, cex.title=0.5, lwd=0.1, collapseTranscripts="meta", showID=FALSE, stacking = "dense")
tetrack<-AnnotationTrack(tes, name = "TEs", stacking = "dense", fill = "lightgreen", shape="box", rotation.title=0, cex.title=0.5, lwd=0.1)

if ( length(htcol) > 0 ) {
	httrack <- HighlightTrack(trackList=tracklist, start = htstart, width = htwidth, col=htcol, fill=htcol2)
	pdf(paste0("plots/",id,"_browser.pdf"),paper="a4")
	plotTracks(list(axistrack, genetrack, tetrack, httrack), sizes=tracksize, title.width=2.5)
	dev.off()
} else {
	tracks<-append(list(axistrack, genetrack, tetrack), tracklist)
	pdf(paste0("plots/",id,"_browser.pdf"),paper="a4")
	plotTracks(tracks, sizes=tracksize, title.width=2.5)
	dev.off()
}





