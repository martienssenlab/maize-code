#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=12G
#$ -l tmp_free=100G
#$ -o ChIPanalysis.log
#$ -j y
#$ -N ChIPanalysis

usage="
##### Script for Maize code ChIP data analysis
#####
##### sh MaiCode_ChIP_analysis.sh -f samplefile [-idr] [-peaks] [-merged] [-h]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Mark, Rep (Rep1, Rep2 or merged), PE or SE
##### 	-h: help, returns usage
##### 
##### It calls broad (for H3K4me1) or narrow (for H3K4me3 and H3K27ac) peaks with Macs2,
##### creates bigwig files (log2 FC vs Input), makes fingerprint plots,
##### makes idr analysis (if both replicates are present in the samplefile), and calculates peak stats
#####
##### If the replicates are to be merged, the Rep value in the samplefile should be 'merged'
#####
##### Requirements: samtools, bedtools, deeptools, macs2, idr
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS

if [ $# -eq 0 ]; then
	printf "$usage\n"
	exit 1
fi

while getopts ":f:sh" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $samplefile ]; then
	printf "Samplefile missing!\n"
	printf "$usage\n"
	exit 1
fi

while read line tissue mark rep paired
do
	#### To merge bam files of replicates if chosen and not already exisiting
	if [[ $rep == "merged" ]]; then
		name=${line}_${tissue}_${mark}_${rep}
		input=${line}_${tissue}_Input_${rep}
		tmpname=${line}_${tissue}_${mark}
		tmpinput=${line}_${tissue}_Input
		if [ ! -s mapped/rmdup_${name}.bam ]; then
			printf "\nMerging replicates of $tmpname\n"
			samtools merge -f -@ $threads mapped/temp_${tmpname}.bam mapped/rmdup_${tmpname}_Rep1.bam mapped/rmdup_${tmpname}_Rep2.bam
			samtools sort -@ $threads -o mapped/rmdup_${tmpname}_merged.bam mapped/temp_${tmpname}.bam
			rm -f mapped/temp_${tmpname}.bam
			samtools index -@ $threads mapped/rmdup_${tmpname}_merged.bam
		fi
		if [ ! -s mapped/rmdup_${input}.bam ]; then
			printf "\nMerging replicates of $tmpinput\n"
			samtools merge -@ $threads mapped/temp_${tmpinput}.bam mapped/rmdup_${tmpinput}_Rep1.bam mapped/rmdup_${tmpinput}_Rep2.bam
			samtools sort -@ $threads -o mapped/rmdup_${tmpinput}_merged.bam mapped/temp_${tmpinput}.bam
			rm -f mapped/temp_${tmpinput}.bam
			samtools index -@ $threads mapped/rmdup_${tmpinput}_merged.bam
		fi
	else
		name=${line}_${tissue}_${mark}_${rep}
		input=${line}_${tissue}_Input_${rep}
	fi
	#### To call either broad or narrow peaks if not already exisiting
	case "$mark" in
		H3K4me1) peaktype="broad";;
		H3K4me3) peaktype="narrow";;
		H3K27ac) peaktype="narrow";;
	esac
	if [[ $paired == "PE" ]]; then
		if [[ $peaktype == "broad" ]] && [ ! -s peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling broad peaks for PE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAMPE -g 2.2e9 -n ${name} --keep-dup "all" --outdir peaks/ --tempdir $TMPDIR --broad
		elif [[ $peaktype == "narrow" ]] && [ ! -s peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling narrow peaks for PE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAMPE -g 2.2e9 -n ${name} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR
		elif [ -s peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\n$peaktype peaks already called for $name\n"
		else
			printf "\nSomething is wrong! Check usage by running the script without arguments\n"
			exit 1
		fi
	elif [[ $paired == "SE" ]]; then
		if [[ $peaktype == "broad" ]] && [ ! -s peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling broad peaks for SE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAM -g 2.2e9 -n ${name} --keep-dup "all" --outdir peaks/ --tempdir $TMPDIR --nomodel --broad
		elif [[ $peaktype == "narrow" ]] && [ ! -s peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling narrow peaks for SE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAM -g 2.2e9 -n ${name} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR --nomodel
		elif [ -s peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\n$peaktype peaks already called for $name\n"
		else
			printf "\nSomething is wrong! Check usage:\n"
			printf "$usage\n"
			exit 1
		fi
	else
		printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
		exit 1
	fi
	#### To create bw files if not already exisiting
	if [ ! -s tracks/${name}.bw ]; then
		printf "\nMaking bigwig files for $name with deeptools version:\n"
		deeptools --version
		bamCompare -b1 mapped/rmdup_${name}.bam -b2 mapped/rmdup_${input}.bam -o tracks/${name}.bw -p $threads --binSize 1 --scaleFactorsMethod "None" --normalizeUsing CPM
	else
		printf "\nBigwig file for $name already exists\n"
	fi
	#### To create fingerprint plots if not already exisiting
	if [ ! -s plots/Fingerprint_${name}.png ]; then
		printf "\nPlotting Fingerprint for $name with deeptools version:\n"
		deeptools --version
		plotFingerprint -b mapped/rmdup_${name}.bam mapped/rmdup_${input}.bam -o plots/Fingerprint_${name}.png -p $threads -l ${name} ${input}
	else
		printf "\nBigwig file for $name already exists\n"
	fi
	touch chkpts/analysis_${name}
done < $samplefile


printf "\nScript finished successfully!\n"

