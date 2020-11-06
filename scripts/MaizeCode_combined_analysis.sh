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
##### Script for Maize code combined data analysis, used by script MaizeCode_analyis.sh if -s was not set and region file exists
#####
##### sh MaiCode_combined_analysis.sh -f samplefile -r regionfile [-s]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Sample, Rep (Rep1, Rep2 or merged), PE or SE
##### 	-r: bedfile containing the regions that are to be ploted over.
##### 	-h: help, returns usage
##### 
##### It produces an Upset plot of the intersection between all ChIP samples, highlighting peaks in the input regions
##### It creates different heatmaps of all ChIP and RNA samples in the samplefile on all input regions
##### Under development:
##### Makes differential peak calling between all pairs of ChIP samples for each mark
##### Call differential expressed genes between all pairs of RNA samples
#####
##### Requirements: bedtools, deeptools, macs2, R (+R packages: ggplot2,readr,UpSetR)
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
while read line tissue sample paired
do
	name=${line}_${tissue}_${sample}
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
		*) datatype="unknown";;
	esac
	bw_list+=("$datatype/tracks/${name}_merged.bw")
	sample_list+=("$name")
	if [[ "$datatype" == "ChIP" ]]; then
		chip_sample_list+=("${name}")
		chip_line_list+=("${line}")
		chip_tissue_list+=("${tissue}")
		chip_mark_list+=("${sample}")
		printf "$line\t$tissue\t$sample\t$paired\n" >> combined/temp_${samplename}_ChIP.txt
	elif [[ "$datatype" == "RNA" ]]; then
		rna_sample_list+=("${name}")
		printf "$line\t$tissue\t$sample\t$paired\n" >> combined/temp_${samplename}_RNA.txt
	else
		printf "\nType of data unknown for ${name}\nSample not processed!\n"
	fi
done < $samplefile

#############################################################################################
########################################### PART2 ###########################################
############################### Overlapping peaks - Upset plot  #############################
#############################################################################################

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
	awk -v OFS="\t" -v s=$sample '{print $1,$2,$3,s}' ChIP/peaks/selected_peaks_${sample}.${peaktype}Peak | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_${samplename}.bed
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
####################### Calling differential peaks accross tissues ##########################
#############################################################################################

# # #### To call differential peaks for all the potential pairs of sample of the same mark between different tissues of the same reference

# uniq_chip_line_list=($(printf "%s\n" "${chip_line_list[@]}" | sort -u))
# uniq_chip_mark_list=($(printf "%s\n" "${chip_mark_list[@]}" | sort -u))

# for line in ${uniq_chip_line_list[@]}
# do
	# sample_line=$(grep "$line" ${chip_sample_list[@]})
	# for mark in ${uniq_chip_mark_list[@]}
	# do
		# sample_line_mark=$(grep "$mark" ${sample_line[@]})
		# numsample=${sample_line_mark[@]}
		# numsamplemin=$((numsample-1))
		# case "$mark" in
			# H3K4me1) 	peaktype="broad"
						# l=500
						# g=400;;
			# H3K4me3) 	peaktype="narrow"
						# l=250
						# g=150;;
			# H3K27ac) 	peaktype="narrow"
						# l=250
						# g=150;;
		# esac
		# if [ $numsample -ge 2 ]; then
			# i=0
			# j=1
			# while [ $i -lt $numsamplemin ]
			# do
				# a1=$(grep "properly paired" reports/flagstat_${sample_line_mark[i]}_Rep1.txt | awk '{print $1}')
				# a2=$(grep "properly paired" reports/flagstat_${sample_line_mark[i]}_Rep2.txt | awk '{print $1}')
				# a=$((a1+a2))
				# while [ $j -le $numsamplemin ]
				# do	
					# b1=$(grep "properly paired" reports/flagstat_${sample_line_mark[j]}_Rep1.txt | awk '{print $1}')
					# b2=$(grep "properly paired" reports/flagstat_${sample_line_mark[j]}_Rep2.txt | awk '{print $1}')
					# b=$((b1+b2))
					# macs2 bdgdiff --t1 ChIP/peaks/${sample_line_mark[i]}_merged_treat_pileup.bdg --c1 ChIP/peaks/${sample_line_mark[i]}_merged_control_lambda.bdg --t2 ChIP/peaks/${sample_line_mark[j]}_merged_treat_pileup.bdg --c2 ChIP/peaks/${sample_line_mark[j]}_merged_control_lambda.bdg --d1 $a --d2 $b -g $g -l $l --outdir combined/peaks --o-prefix diff_${sample_line_mark[i]}_${sample_line_mark[j]}
					# j=$((j+1))
				# done
				# i=$((i+1))
			# done
		# fi
	# done
# done


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
printf "\nComputing scale-regions matrix for $analysisname\n"
computeMatrix scale-regions --missingDataAsZero --skipZeros -R $regionfile -S ${bw_list[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/regions_${analysisname}.gz
printf "\nComputing reference-point on TSS matrix for $analysisname\n"
computeMatrix reference-point --referencePoint "TSS" --missingDataAsZero --skipZeros -R $regionfile -S ${bw_list[@]} -bs 50 -b 2000 -a 8000 -p $threads -o combined/matrix/tss_${analysisname}.gz
printf "\nComputing reference-point on TES matrix for $analysisname\n"
computeMatrix reference-point --referencePoint "TES" --missingDataAsZero --skipZeros -R $regionfile -S ${bw_list[@]} -bs 50 -b 8000 -a 2000 -p $threads -o combined/matrix/tes_${analysisname}.gz

for matrix in regions tss tes
do
	printf "\nGetting scales (10th and 90th quantiles) for ${matrix} matrix\n"
	computeMatrixOperations dataRange -m combined/matrix/${matrix}_${analysisname}.gz > combined/matrix/values_${matrix}_${analysisname}.txt
	mins=()
	maxs=()
	for sample in ${sample_list[@]}
	do
		mini=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $5}')
		mins+=("$mini")
		maxi=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $6}')
		maxs+=("$maxi")
	done
	printf "\nPlotting heatmap for by $matrix of $analysisname\n"
	plotHeatmap -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_${matrix}.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sample_list[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]}
	printf "\nPlotting heatmap for by $matrix of $analysisname in 5 clusters (kmeans)\n"
	plotHeatmap -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_${matrix}_k5.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sample_list[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --kmeans 5 --outFileSortedRegions combined/matrix/${analysisname}_${matrix}_regions_k5.txt
done

printf "\nCombined analysis script finished successfully\n"
touch combined/chkpts/${analysisname}
