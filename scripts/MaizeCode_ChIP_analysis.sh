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
##### sh MaiCode_ChIP_analysis.sh -f samplefile -r regionfile [-s]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Mark, Rep (Rep1, Rep2 or merged), PE or SE
##### 	-r: regionfile (bed file) containing the regions that want to be ploted over
##### 		(safest to use a full path to the region file)
##### 	-h: help, returns usage
##### 
##### It calls broad (for H3K4me2) or narrow (for H3K4me3 and H3K27ac) peaks with Macs2,
##### creates bigwig files (log2 FC vs Input), makes fingerprint plots,
##### idr analysis (if both replciates are present in samplefile), calculates peak stats,
##### and creates an upset plots to see overlap between samples, highlighting peaks in gene bodies
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

while getopts ":f:r:sh" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		r)	export regionfile=${OPTARG};;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $samplefile ] || [ ! $regionfile ]; then
	printf "Missing arguments!\n"
	printf "$usage\n"
	exit 1
fi

if [ ! -d ./peaks ]; then
	mkdir ./peaks
fi

if [ ! -d ./plots ]; then
	mkdir ./plots
fi

line_list=()
tissue_list=()
mark_list=()
sample_list=()
bw_list=()
bam_list=()
while read line tissue mark rep paired
do
	#### To merge bam files of replicates if chosen and not already exisiting
	if [[ $rep == "merged" ]]; then
		name=${line}_${tissue}_${mark}_${rep}
		input=${line}_${tissue}_Input_${rep}
		tmpname=${line}_${tissue}_${mark}
		tmpinput=${line}_${tissue}_Input
		if [ ! -f mapped/rmdup_${name}.bam ]; then
			printf "\nMerging replicates of $tmpname\n"
			samtools merge -f -@ $threads mapped/temp_${tmpname}.bam mapped/rmdup_${tmpname}_Rep1.bam mapped/rmdup_${tmpname}_Rep2.bam
			samtools sort -@ $threads -o mapped/rmdup_${tmpname}_merged.bam mapped/temp_${tmpname}.bam
			rm -f mapped/temp_${tmpname}.bam
			samtools index -@ $threads mapped/rmdup_${tmpname}_merged.bam
		fi
		if [ ! -f mapped/rmdup_${input}.bam ]; then
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
	#### To append the lists (used for labels and other stuff)
	sample_list+=("$name")
	line_list+=("$line")
	tissue_list+=("$tissue")
	mark_list+=("$mark")
	bam_list+=("mapped/rmdup_${name}.bam" "mapped/rmdup_${input}.bam")	
	#### To call either broad or narrow peaks if not already exisiting
	case "$mark" in
		H3K4me1) peaktype="broad";;
		H3K4me3) peaktype="narrow";;
		H3K27ac) peaktype="narrow";;
	esac
	if [[ $paired == "PE" ]]; then
		if [[ $peaktype == "broad" ]] && [ ! -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling broad peaks for PE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAMPE -g 2.2e9 -n ${name} --keep-dup "all" --outdir peaks/ --tempdir $TMPDIR --broad
		elif [[ $peaktype == "narrow" ]] && [ ! -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling narrow peaks for PE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAMPE -g 2.2e9 -n ${name} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR
		elif [ -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\n$peaktype peaks already called for $name\n"
		else
			printf "\nSomething is wrong! Check usage by running the script without arguments\n"
			exit 1
		fi
	elif [[ $paired == "SE" ]]; then
		if [[ $peaktype == "broad" ]] && [ ! -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling broad peaks for SE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAM -g 2.2e9 -n ${name} --keep-dup "all" --outdir peaks/ --tempdir $TMPDIR --nomodel --broad
		elif [[ $peaktype == "narrow" ]] && [ ! -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling narrow peaks for SE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAM -g 2.2e9 -n ${name} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR --nomodel
		elif [ -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\n$peaktype peaks already called for $name\n"
		else
			printf "\nSomething is wrong! Check usage by running the script without arguments\n"
			exit 1
		fi
	else
		printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
		exit 1
	fi
	#### To create bw files if not already exisiting
	if [ ! -f tracks/${name}.bw ]; then
		printf "\nMaking bigwig files for $name with deeptools version:\n"
		deeptools --version
		bamCompare -b1 mapped/rmdup_${name}.bam -b2 mapped/rmdup_${input}.bam -o tracks/${name}.bw -p $threads --binSize 1 --scaleFactorsMethod "None" --normalizeUsing CPM
	else
		printf "\nBigwig file for $name already exists\n"
	fi
	#### To create fingerprint plots if not already exisiting
	if [ ! -f plots/Fingerprint_${name}.png ]; then
		printf "\nPlotting Fingerprint for $name with deeptools version:\n"
		deeptools --version
		plotFingerprint -b mapped/rmdup_${name}.bam mapped/rmdup_${input}.bam -o plots/Fingerprint_${name}.png -p $threads -l ${name} ${input}
	else
		printf "\nBigwig file for $name already exists\n"
	fi
	#### To append the list of bw files to be used by deeptools
	bw_list+=("tracks/${name}.bw")
done < $samplefile

if [[ $keepgoing == "STOP" ]]; then
	printf "\nScript finished successfully without analysis\n"		
	exit 0
fi	

tmp1=${samplefile##*temp_}
samplename=${tmp1%_ChIP*}

#### To get the list of unique information (i.e. get rid of potential repeats)

uniq_line_list=($(printf "%s\n" "${line_list[@]}" | sort -u))
uniq_tissue_list=($(printf "%s\n" "${tissue_list[@]}" | sort -u))
uniq_mark_list=($(printf "%s\n" "${mark_list[@]}" | sort -u))

#### To create idr analysis and summary peak statistics

printf "Line\tTissue\Mark\tPeaks_in_rep1\tPeaks_in_Rep2\tCommon_peaks\tCommon_peaks_IDR<=0.05\n" > reports/summary_peaks_${samplename}.txt
for line in ${uniq_line_list[@]}
do
	for tissue in ${uniq_tissue_list[@]}
	do
		for mark in ${uniq_mark_list[@]}
		do
			if [[ " ${sample_list[@]} " =~ "${line}_${tissue}_${mark}" ]]; then
				case "$mark" in
					H3K4me1) peaktype="broad";;
					H3K4me3) peaktype="narrow";;
					H3K27ac) peaktype="narrow";;
				esac
				#### To get IDR analysis on replicates (if they were both present in the samplefile)
				if [[ " ${sample_list[@]} " =~ " ${line}_${tissue}_${mark}_Rep1 " ]] && [[ " ${sample_list[@]} " =~ " ${line}_${tissue}_${mark}_Rep2 " ]]; then
					if [ ! -f peaks/idr_${line}_${tissue}_${mark}.${peaktype}Peak ]; then
						printf "\nDoing IDR analysis on both replicates from ${line}_${tissue}_${mark} ($peaktype peaks) with idr version:\n"
						idr --version
						idr --input-file-type ${peaktype}Peak --output-file-type ${peaktype}Peak --samples peaks/${line}_${tissue}_${mark}_Rep1_peaks.${peaktype}Peak peaks/${line}_${tissue}_${mark}_Rep2_peaks.${peaktype}Peak -o peaks/idr_${line}_${tissue}_${mark}.${peaktype}Peak -l reports/idr_${line}_${tissue}_${mark}.log --plot
					else
						printf "\nIDR analysis already done for ${line}_${tissue}_${mark}\n"
					fi
					#### To get some peaks stats for each type of mark (also if both replicates were present)
					printf "\nCalculating peak stats for ${line}_${tissue}_${mark} in ${peaktype} peaks\n"
					rep1=$(awk '{print $1,$2,$3}' peaks/${line}_${tissue}_${mark}_Rep1_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
					rep2=$(awk '{print $1,$2,$3}' peaks/${line}_${tissue}_${mark}_Rep2_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
					common=$(awk '{print $1,$2,$3}' peaks/idr_${line}_${tissue}_${mark}.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
					idr=$(awk '$5>=540 {print $1,$2,$3}' peaks/idr_${line}_${tissue}_${mark}.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
					awk -v OFS="\t" -v a=$line -v b=$tissue -v c=$mark -v d=$rep1 -v e=$rep2 -v f=$common -v g=$idr 'BEGIN {print a,b,c,d,e,f" ("f/d*100"%rep1;"f/e*100"%rep2)",g" ("g/f*100"%common)"}' >> reports/summary_peaks_${samplename}.txt
				fi
			fi
		done
	done
done

#### To make a single file containing all overlapping peaks

printf "\nPreparing merged peaks file for $samplename\n"
if [ -f peaks/tmp_peaks_${samplename}.bed ]; then
	rm -f peaks/tmp_peaks_${samplename}.bed
fi

for sample in ${sample_list[@]}
do
	case "$sample" in
		*H3K4me1*) peaktype="broad";;
		*H3K4me3*) peaktype="narrow";;
		*H3K27ac*) peaktype="narrow";;
	esac
	awk -v OFS="\t" -v s=$sample '{print $1,$2,$3,s}' peaks/${sample}_peaks.${peaktype}Peak | uniq >> peaks/tmp_peaks_${samplename}.bed
done
sort -k1,1 -k2,2n peaks/tmp_peaks_${samplename}.bed > peaks/tmp2_peaks_${samplename}.bed
bedtools merge -i peaks/tmp2_peaks_${samplename}.bed -c 4 -o distinct | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2,$3,"Peak_"NR,$4}'> peaks/tmp3_peaks_${samplename}.bed

#### To get distance to closest gene (and the gene model name)

printf "\nGetting closest region of $samplename to $regionfile\n"
bedtools closest -a peaks/tmp3_peaks_${samplename}.bed -b $regionfile -D ref | awk -v OFS="\t" '{print $1,$2,$3,$4,$12,".",$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$9,$7}' > peaks/peaks_${samplename}.bed
rm -f peaks/tmp*

#### To create a matrix of peak presence in each sample

printf "\nCreating matrix file for $samplename\n"
for sample in ${sample_list[@]}
do
	printf "$sample\n" > peaks/temp_col_${sample}.txt
	awk -v OFS="\t" -v s=$sample '{if ($0 ~ s) print "1"; else print "0"}' peaks/peaks_${samplename}.bed >> peaks/temp_col_${sample}.txt
done

#### To group peaks based on their distance (gene body (x=0), promoter (0<x<2kb upstream), terminator (0<x<2kb downstream), distal)

awk -v OFS="\t" 'BEGIN {printf "PeakID\tDistance\n"} {if ($5<-2000) d="Distal"; else if ($5<0) d="Promoter"; else if ($5==0) d="Gene_body"; else if ($5>2000) d="Distal"; else d="Terminator"; print $4,d}' peaks/peaks_${samplename}.bed > peaks/temp_col_AAA.txt

paste peaks/temp_col_*.txt | uniq > peaks/matrix_upset_${samplename}.txt

rm peaks/temp_col_*.txt

#### To make an Upset plot highlighting peaks in gene bodies

printf "\nCreating Upset plot for $samplename with R version:\n"
R --version
Rscript --vanilla ~/data/Scripts/MaizeCode_R_Upset.r peaks/matrix_upset_${samplename}.txt ${samplename}

printf "\nScript finished successfully!\n"		
