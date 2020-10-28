#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=12G
#$ -l tmp_free=200G
#$ -o combinedanalysis.log
#$ -j y
#$ -N combinedanalysis

usage="
##### Script for Maize code combined data analysis
#####
##### sh MaiCode_analysis.sh -f samplefile -r regionfile [-s]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Sample, Rep (Rep1, Rep2 or merged), PE or SE
##### 	-r: bedfiles containing the regions that want to be ploted over
##### 		(safest to use a full path to the region file)
##### 	-h: help, returns usage
##### 
##### It sends each type of sample to its specific analysis file (MaizeCode_ChIP_analysis.sh or MaizeCode_RNA_analysis.sh)
##### Then combines into deeper analysis (to detail...)
#####
##### Requirements: samtools, bedtools, deeptools, macs2, idr, R (+R packages: ggplot2,readr,UpSetR)
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS
# # export mc_dir=$(dirname "$0")
export mc_dir="${HOME}/data/Scripts/MaizeCode/"
printf "\nRunning MaizeCode scripts from ${mc_dir} in working directory ${PWD}\n"

if [ $# -eq 0 ]; then
	printf "$usage\n"
	exit 1
fi

while getopts ":f:r:h" opt; do
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

if [ ! $samplefile ]; then
	printf "Missing samplefile!\n"
	printf "$usage\n"
	exit 1
fi

if [ ! $regionfile ]; then
	printf "Missing regionfile\n"
	printf "$usage\n"
	exit 1
fi

#############################################################################################
########################################### PART1 ###########################################
############################## Prepare the ChIP and RNA samples #############################
#############################################################################################

tmp1=${samplefile##*/}
samplename=${tmp1%_analysis*}
tmp2=${regionfile##*/}
regionname=${tmp2%.*}
analysisname=${samplename}_on_${regionname}

printf "\nStarting analysis: $analysisname\n"

if [ -s combined/temp_${samplename}_ChIP.txt ]; then
	rm -f combined/temp_${samplename}_ChIP.txt
fi

if [ -s combined/temp_${samplename}_RNA.txt ]; then
	rm -f combined/temp_${samplename}_RNA.txt
fi

#### To append the lists (used for labels and other stuff)

sample_list=()
chip_sample_list=()
chip_line_list=()
chip_tissue_list=()
chip_mark_list=()
rna_sample_list=()
bw_list=()
while read line tissue sample rep paired
do
	name=${line}_${tissue}_${sample}_${rep}
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
		*) datatype="unknown";;
	esac
	bw_list+=("$datatype/tracks/${name}.bw")
	sample_list+=("$name")
	if [[ "$datatype" == "ChIP" ]]; then
		chip_sample_list+=("${name}")
		chip_line_list+=("${line}")
		chip_tissue_list+=("${tissue}")
		chip_mark_list+=("${sample}")
		printf "$line\t$tissue\t$sample\t$rep\t$paired\n" >> combined/temp_${samplename}_ChIP.txt
	elif [[ "$datatype" == "RNA" ]]; then
		rna_sample_list+=("${name}")
		printf "$line\t$tissue\t$sample\t$rep\t$paired\n" >> combined/temp_${samplename}_RNA.txt
	else
		printf "\nType of data unknown for ${name}\nSample not processed!\n"
	fi
done < $samplefile

#############################################################################################
########################################### PART2 ###########################################
################################## Combine the ChIP samples #################################
#############################################################################################

#### To create idr analysis and summary peak statistics

uniq_chip_line_list=($(printf "%s\n" "${chip_line_list[@]}" | sort -u))
uniq_chip_tissue_list=($(printf "%s\n" "${chip_tissue_list[@]}" | sort -u))
uniq_chip_mark_list=($(printf "%s\n" "${chip_mark_list[@]}" | sort -u))

printf "Line\tTissue\tMark\tPeaks_in_rep1\tPeaks_in_Rep2\tCommon_peaks\tCommon_peaks_IDR<=0.05\n" > combined/peaks/summary_peaks_${samplename}.txt

for line in ${uniq_chip_line_list[@]}
do
	for tissue in ${uniq_chip_tissue_list[@]}
	do
		for mark in ${uniq_chip_mark_list[@]}
		do
			if [[ " ${chip_sample_list[@]} " =~ "${line}_${tissue}_${mark}" ]]; then
				case "$mark" in
					H3K4me1) peaktype="broad";;
					H3K4me3) peaktype="narrow";;
					H3K27ac) peaktype="narrow";;
				esac
				#### To get IDR analysis on replicates (if they were both present in the samplefile)
				if [[ " ${chip_sample_list[@]} " =~ " ${line}_${tissue}_${mark}_Rep1 " ]] && [[ " ${chip_sample_list[@]} " =~ " ${line}_${tissue}_${mark}_Rep2 " ]]; then
					if [ ! -s ChIP/peaks/idr_${line}_${tissue}_${mark}.${peaktype}Peak ]; then
						printf "\nDoing IDR analysis on both replicates from ${line}_${tissue}_${mark} ($peaktype peaks) with idr version:\n"
						idr --version
						idr --input-file-type ${peaktype}Peak --output-file-type ${peaktype}Peak --samples ChIP/peaks/${line}_${tissue}_${mark}_Rep1_peaks.${peaktype}Peak ChIP/peaks/${line}_${tissue}_${mark}_Rep2_peaks.${peaktype}Peak -o ChIP/peaks/idr_${line}_${tissue}_${mark}.${peaktype}Peak -l ChIP/reports/idr_${line}_${tissue}_${mark}.log --plot
					else
						printf "\nIDR analysis already done for ${line}_${tissue}_${mark}\n"
					fi
					#### To get some peaks stats for each type of mark (also if both replicates were present)
					printf "\nCalculating peak stats for ${line}_${tissue}_${mark} in ${peaktype} peaks\n"
					rep1=$(awk '{print $1,$2,$3}' ChIP/peaks/${line}_${tissue}_${mark}_Rep1_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
					rep2=$(awk '{print $1,$2,$3}' ChIP/peaks/${line}_${tissue}_${mark}_Rep2_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
					common=$(awk '{print $1,$2,$3}' ChIP/peaks/idr_${line}_${tissue}_${mark}.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
					idr=$(awk '$5>=540 {print $1,$2,$3}' ChIP/peaks/idr_${line}_${tissue}_${mark}.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
					awk -v OFS="\t" -v a=$line -v b=$tissue -v c=$mark -v d=$rep1 -v e=$rep2 -v f=$common -v g=$idr 'BEGIN {print a,b,c,d,e,f" ("f/d*100"%rep1;"f/e*100"%rep2)",g" ("g/f*100"%common)"}' >> combined/peaks/summary_peaks_${samplename}.txt
				fi
			fi
		done
	done
done

#### To make a single file containing all overlapping peaks

printf "\nPreparing merged peaks file for $samplename\n"
if [ -s combined/peaks/tmp_peaks_${samplename}.bed ]; then
	rm -f combined/peaks/tmp_peaks_${samplename}.bed
fi

for sample in ${chip_sample_list[@]}
do
	case "$sample" in
		*H3K4me1*) peaktype="broad";;
		*H3K4me3*) peaktype="narrow";;
		*H3K27ac*) peaktype="narrow";;
	esac
	awk -v OFS="\t" -v s=$sample '{print $1,$2,$3,s}' ChIP/peaks/${sample}_peaks.${peaktype}Peak | uniq >> combined/peaks/tmp_peaks_${samplename}.bed
done
sort -k1,1 -k2,2n combined/peaks/tmp_peaks_${samplename}.bed > combined/peaks/tmp2_peaks_${samplename}.bed
bedtools merge -i combined/peaks/tmp2_peaks_${samplename}.bed -c 4 -o distinct | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2,$3,"Peak_"NR,$4}'> combined/peaks/tmp3_peaks_${samplename}.bed

#### To get distance to closest gene (and the gene model name)

printf "\nGetting closest region of $samplename to $regionfile\n"
bedtools closest -a combined/peaks/tmp3_peaks_${samplename}.bed -b $regionfile -D ref | awk -v OFS="\t" '{print $1,$2,$3,$4,$12,".",$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$9,$7}' > combined/peaks/peaks_${analysisname}.bed
rm -f combined/peaks/tmp*

#### To create a matrix of peak presence in each sample

printf "\nCreating matrix file for $samplename\n"
for sample in ${chip_sample_list[@]}
do
	printf "$sample\n" > combined/peaks/temp_col_${sample}.txt
	awk -v OFS="\t" -v s=$sample '{if ($0 ~ s) print "1"; else print "0"}' combined/peaks/peaks_${analysisname}.bed >> combined/peaks/temp_col_${sample}.txt
done

#### To group peaks based on their distance (gene body (x=0), promoter (0<x<2kb upstream), terminator (0<x<2kb downstream), distal)

awk -v OFS="\t" 'BEGIN {printf "PeakID\tDistance\n"} {if ($5<-2000) d="Distal"; else if ($5<0) d="Promoter"; else if ($5==0) d="Gene_body"; else if ($5>2000) d="Distal"; else d="Terminator"; print $4,d}' combined/peaks/peaks_${analysisname}.bed > combined/peaks/temp_col_AAA.txt

paste combined/peaks/temp_col_*.txt | uniq > combined/peaks/matrix_upset_${samplename}.txt

rm -f combined/peaks/temp_col_*.txt

#### To make an Upset plot highlighting peaks in gene bodies

printf "\nCreating Upset plot for $samplename with R version:\n"
R --version
Rscript --vanilla ${mc_dir}/MaizeCode_R_Upset.r combined/peaks/matrix_upset_${samplename}.txt ${analysisname}


#############################################################################################
########################################### PART3 ###########################################
################################## Combine the RNA samples ##################################
#############################################################################################




#############################################################################################
########################################### PART4 ###########################################
####################################### Making heatmaps  ####################################
#############################################################################################

#### To make heatmaps and profiles with deeptools
#### By default, it does both scale-regions and reference-point on start of bedfile provided
#### By default, it does heatmap on all the data, heatmap with 5 kmeans, and corresponding profiles
#### Probably need to edit many parameters depending on the purpose of the analysis

printf "\nDoing analysis for $analysisname with deeptools version:\n"
deeptools --version

#### Computing the matrix
if [ ! -f combined/matrix/regions_${analysisname}.gz ]; then
	printf "\nComputing scale-regions matrix for $analysisname\n"
	computeMatrix scale-regions -R $regionfile -S ${bw_list[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/regions_${analysisname}.gz
fi
if [ ! -f combined/matrix/tss_${analysisname}.gz ]; then
	printf "\nComputing reference-point on TSS matrix for $analysisname\n"
	computeMatrix reference-point --referencePoint "TSS" -R $regionfile -S ${bw_list[@]} -bs 50 -b 2000 -a 6000 -p $threads -o combined/matrix/tss_${analysisname}.gz
fi

#### Ploting heatmaps
printf "\nPlotting full heatmap for scale-regions of $analysisname\n"
plotHeatmap -m combined/matrix/regions_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_regions.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sample_list[@]} --colorMap 'seismic'
printf "\nPlotting heatmap for scale-regions of $analysisname split in 3 kmeans\n"
plotHeatmap -m combined/matrix/regions_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_regions_k3.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sample_list[@]} --colorMap 'seismic' --kmeans 3 --outFileSortedRegions combined/matrix/${analysisname}_sortedregions_k3.txt

printf "\nPlotting full heatmap for reference-point TSS of $analysisname\n"
plotHeatmap -m combined/matrix/tss_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_tss.pdf --sortRegions descend --sortUsing region_length --samplesLabel ${sample_list[@]} --colorMap 'seismic'
printf "\nPlotting heatmap for reference-point TSS of $analysisname split in 3 kmeans\n"
plotHeatmap -m combined/matrix/tss_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_tss_k3.pdf --sortRegions descend --sortUsing region_length --samplesLabel ${sample_list[@]} --colorMap 'seismic' --kmeans 3 --outFileSortedRegions combined/matrix/${analysisname}_sortedtss_k3.txt

# #### Plotting Metaplot profiles
# printf "\nPlotting metaplot profiles for scale-regions of $analysisname\n"
# plotProfile -m combined/matrix/regions_${analysisname}.gz -out combined/plots/${analysisname}_profiles_regions.pdf --plotType lines --averageType mean --perGroup --samplesLabel ${sample_list[@]}
# printf "\nPlotting metaplot profiles for scale-regions of $analysisname split in 5 kmeans\n"
# plotProfile -m combined/matrix/regions_${analysisname}.gz -out combined/plots/${analysisname}_profiles_regions_k5.pdf --plotType lines --averageType mean --perGroup --samplesLabel ${sample_list[@]} --kmeans 5

# printf "\nPlotting metaplot profiles for reference-point TSS of $analysisname\n"
# plotProfile -m combined/matrix/tss_${analysisname}.gz -out combined/plots/${analysisname}_profiles_tss.pdf --plotType lines --averageType mean --perGroup --samplesLabel ${sample_list[@]}
# printf "\nPlotting metaplot profiles for reference-point TSS of $analysisname split in 5 kmeans\n"
# plotProfile -m combined/matrix/tss_${analysisname}.gz -out combined/plots/${analysisname}_profiles_tss_k5.pdf --plotType lines --averageType mean --perGroup --samplesLabel ${sample_list[@]} --kmeans 5

printf "\nCombined analysis script finished successfully\n"
touch combined/chkpts/${analysisname}
