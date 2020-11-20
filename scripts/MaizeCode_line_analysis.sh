#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=12G
#$ -l tmp_free=200G
#$ -o lineanalysis.log
#$ -j y
#$ -N lineanalysis

usage="
##### Script for Maize code inbred line data analysis, used by script MaizeCode_analyis.sh if -s was not set and region file exists
#####
##### sh MaiCode_line_analysis.sh -f samplefile -r regionfile [-s]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Sample, PE or SE, Reference genome directory
##### 	-r: bedfile containing the regions that are to be ploted over
##### 	-h: help, returns usage
##### 
##### It produces an Upset plot of the intersection between all ChIP samples, highlighting peaks in the input regions
##### It creates different heatmaps of all ChIP and RNA samples in the samplefile on all input regions
##### Under development:
##### Differential peak calling between all pairs of ChIP samples for each mark
##### Differential gene expression analysis between all pairs of RNA samples
##### Differential peak/TSS calling between all pairs of RAMPAGE samples
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
regionname=${tmp2%%.*}
analysisname=${samplename}_on_${regionname}

printf "\nStarting analysis: $analysisname\n"

#### To append the lists (used for labels and other stuff)

chip_sample_list=()
chip_tissue_list=()
chip_mark_list=()
chip_bw_list=()
rna_type_list=()
rnaseq_sample_list=()
rnaseq_tissue_list=()
rampage_sample_list=()
rampage_tissue_list=()
shrna_sample_list=()
shrna_tissue_list=()
rnaseq_bw_list_plus=()
rnaseq_bw_list_minus=()
while read line tissue sample paired ref_dir
do
	name=${line}_${tissue}_${sample}
	case "$sample" in
		H*|Input) export datatype="ChIP";;
		*RNA*|RAMPAGE) export datatype="RNA";;
		*) export datatype="unknown";;
	esac
	ref=${ref_dir##*/}
	if [[ ! "${ref_list[@]}" =~ "${ref}" ]]; then
		ref_list+=("$ref")
	fi
	if [[ "$datatype" == "ChIP" ]]; then
		chip_bw_list+=("$datatype/tracks/${name}_merged.bw")
		chip_sample_list+=("${name}")
		chip_tissue_list+=("${tissue}")
		chip_mark_list+=("${sample}")
	elif [[ "$datatype" == "RNA" ]]; then
		rna_type_list+=("${sample}")
		if [[ "$sample" == "RNAseq" ]]; then
			rnaseq_bw_list_plus+=("$datatype/tracks/${name}_merged_plus.bw")
			rnaseq_bw_list_minus+=("$datatype/tracks/${name}_merged_minus.bw")
			rnaseq_sample_list+=("${name}")
			rnaseq_tissue_list+=("${tissue}")
		elif [[ "$sample" == "RAMPAGE" ]]; then
			rampage_sample_list+=("${name}")
			rampage_tissue_list+=("${tissue}")
		elif [[ "$sample" == "shRNA" ]]; then
			shrna_sample_list+=("${name}")
			shrna_tissue_list+=("${tissue}")
		else 
			printf "\nType of RNA sample unknown\n"
		fi
	else
		printf "\nType of data unknown for ${name}\nSample not processed!\n"
	fi
done < $samplefile

if [ ! ${#ref_list[@]} -eq 1 ]; then
	printf "\nThere are multiple references in the samplefile! This analysis cannot be performed!\n"
	exit 1
else
	export ref=${ref_list[0]}
fi

#############################################################################################
########################################### PART2 ###########################################
########################## Overlapping ChIPseq peaks - Upset plot  ##########################
#############################################################################################

#### To make a single file containing all overlapping peaks

if [ ${#chip_sample_list[@]} -ge 1 ]; then
	printf "\nPreparing merged peaks file for $analysisname\n"
	if [ -s combined/peaks/tmp_peaks_${analysisname}.bed ]; then
		rm -f combined/peaks/tmp_peaks_${analysisname}.bed
	fi
	for sample in ${chip_sample_list[@]}
	do
		case "$sample" in
			*H3K4me1*) export peaktype="broad";;
			*H3K4me3*) export peaktype="narrow";;
			*H3K27ac*) export peaktype="narrow";;
		esac
		awk -v OFS="\t" -v s=$sample '{print $1,$2,$3,s}' ChIP/peaks/selected_peaks_${sample}.${peaktype}Peak | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_${analysisname}.bed
	done
	sort -k1,1 -k2,2n combined/peaks/tmp_peaks_${analysisname}.bed > combined/peaks/tmp2_peaks_${analysisname}.bed
	bedtools merge -i combined/peaks/tmp2_peaks_${analysisname}.bed -c 4 -o distinct | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2,$3,"Peak_"NR,$4}'> combined/peaks/tmp3_peaks_${analysisname}.bed
	#### To get distance to closest gene (and the gene model name)
	printf "\nGetting closest region of $samplename to $regionfile\n"
	bedtools closest -a combined/peaks/tmp3_peaks_${analysisname}.bed -b $regionfile -D ref | awk -v OFS="\t" '{print $1,$2,$3,$4,$12,".",$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$9,$7}' > combined/peaks/peaks_${analysisname}.bed
	rm -f combined/peaks/tmp*_peaks_${analysisname}.bed
	#### To create a matrix of peak presence in each sample
	printf "\nCreating matrix file for $samplename\n"
	for sample in ${chip_sample_list[@]}
	do
		printf "$sample\n" > combined/peaks/temp_col_${analysisname}_${sample}.txt
		awk -v OFS="\t" -v s=$sample '{if ($0 ~ s) print "1"; else print "0"}' combined/peaks/peaks_${analysisname}.bed >> combined/peaks/temp_col_${analysisname}_${sample}.txt
	done
	#### To group peaks based on their distance (gene body (x=0), promoter (0<x<2kb upstream), terminator (0<x<2kb downstream), distal)
	awk -v OFS="\t" 'BEGIN {printf "PeakID\tDistance\n"} {if ($5<-2000) d="Distal"; else if ($5<0) d="Promoter"; else if ($5==0) d="Gene_body"; else if ($5>2000) d="Distal"; else d="Terminator"; print $4,d}' combined/peaks/peaks_${analysisname}.bed > combined/peaks/temp_col_${analysisname}_AAA.txt
	paste combined/peaks/temp_col_${analysisname}_*.txt | uniq > combined/peaks/matrix_upset_${analysisname}.txt
	rm -f combined/peaks/temp_col_${analysisname}_*.txt
	#### To make an Upset plot highlighting peaks in gene bodies
	printf "\nCreating Upset plot for $samplename with R version:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_Upset.r combined/peaks/matrix_upset_${analysisname}.txt ${analysisname}
fi

############################################################################################
########################################## PART3 ###########################################
############################# Overlapping TSS - Upset plot  ################################
############################################################################################

#### To make a single file containing all overlapping TSS

############################################################################################
########################################## PART4 ###########################################
################# Differential gene expression analysis between tissues ####################
############################################################################################

if [ ${#rnaseq_tissue_list[@]} -ge 2 ]; then
	#### To make a count table for all RNAseq samples in samplefile
	printf "\nPreparing count table for RNAseq samples in $analysisname\n"
	printf "Replicate\tSample\n" > combined/DEG/samples_${analysisname}.txt
	i=0
	for sample in ${rnaseq_sample_list[@]}
	do
		tissuei=${rnaseq_tissue_list[i]}
		printf "${tissuei}_Rep1\t${tissuei}\n${tissuei}_Rep2\t${tissuei}\n" >> combined/DEG/samples_${analysisname}.txt
		grep "gene:" RNA/mapped/map_${sample}_Rep1_ReadsPerGene.out.tab | sed 's/gene://' | awk -v OFS="\t" -v t=$tissuei 'BEGIN {print t"_Rep1"} {print $2}' > combined/DEG/col_AZ_${i}_${analysisname}_${sample}_Rep1.txt
		grep "gene:" RNA/mapped/map_${sample}_Rep2_ReadsPerGene.out.tab | sed 's/gene://' | awk -v OFS="\t" -v t=$tissuei 'BEGIN {print t"_Rep2"} {print $2}' > combined/DEG/col_AZ_${i}_${analysisname}_${sample}_Rep2.txt
		if [ $i -eq 0 ]; then
			grep "gene:" RNA/mapped/map_${sample}_Rep1_ReadsPerGene.out.tab | sed 's/gene://' | awk -v OFS="\t" 'BEGIN {print "gene_ID"} {print $1}' > combined/DEG/col_AA_0_${analysisname}.txt
		fi
		i=$((i+1))
	done
	paste combined/DEG/col_A*_${analysisname}* > combined/DEG/counts_${analysisname}.txt
	rm -f combined/DEG/col_A*_${analysisname}*
	#### To run the DEG analysis on R
	printf "\nLaunching DEG analysis with R version:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_DEG.r combined/DEG/counts_${analysisname}.txt combined/DEG/samples_${analysisname}.txt ${analysisname} $regionfile
	#### To extract DEG only called in one tissue
	if [ ${#rnaseq_tissue_list[@]} -ge 3 ]; then
		for file in $(ls combined/DEG/DEG_${analysisname}_*vs*.txt)
		do
			awk -v OFS="\t" '{print $1,$2,$3,$4,".",$6}' ${file} > ${file}.temp.bed
		done
		for tissue in ${rnaseq_tissue_list[@]}
		do
			printf "\nGetting DEGs specific in $tissue\n"
			filearray=( $(ls combined/DEG/DEG_${analysisname}_*vs*.txt.temp.bed | grep "$tissue") )
			i=0
			max=$((${#filearray[@]}-1))
			cat ${filearray[i]} > combined/DEG/temp_tissue_spec_DEG_${i}.txt
			while [ $i -lt $max ]
			do
				j=$((i+1))
				bedtools intersect -wa -a combined/DEG/temp_tissue_spec_DEG_${i}.txt -b ${filearray[j]} > combined/DEG/temp_tissue_spec_DEG_${j}.txt
				i=$((i+1))
			done
			cat combined/DEG/temp_tissue_spec_DEG_${max}.txt > combined/DEG/only_${tissue}_DEG_${analysisname}.bed
		done
		rm -f combined/DEG/DEG_${analysisname}_*vs*.txt.temp.bed
		rm -f combined/DEG/temp_tissue_spec_DEG_*
	fi
else 
	printf "\nNo differential gene expression analysis performed (not enough samples)\n"
fi


############################################################################################
########################################## PART5 ###########################################
##################### Differential peaks analysis between tissues ##########################
############################################################################################

# #### To call differential peaks for all pairs of tissues of the same mark

# uniq_chip_mark_list=($(printf "%s\n" "${chip_mark_list[@]}" | sort -u))
# for mark in ${uniq_chip_mark_list[@]}
# do
	# sample_line_mark=()
	# for sample in ${chip_sample_list[@]}
	# do
		# if [[ $sample =~ $mark ]]; then
			# sample_line_mark+=("$sample")
		# fi
	# done
	# numsample=${#sample_line_mark[@]}
	# numsamplemin1=$((numsample-1))
	# numsamplemin2=$((numsample-2))
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
		# while [ $i -lt $numsamplemin2 ]
		# do
			# a1=$(grep "properly paired" ChIP/reports/flagstat_${sample_line_mark[i]}_Rep1.txt | awk '{print $1}')
			# a2=$(grep "properly paired" ChIP/reports/flagstat_${sample_line_mark[i]}_Rep2.txt | awk '{print $1}')
			# a=$((a1+a2))
			# j=$((i+1))
			# while [ $j -le $numsamplemin1 ]
			# do	
				# b1=$(grep "properly paired" ChIP/reports/flagstat_${sample_line_mark[j]}_Rep1.txt | awk '{print $1}')
				# b2=$(grep "properly paired" ChIP/reports/flagstat_${sample_line_mark[j]}_Rep2.txt | awk '{print $1}')
				# b=$((b1+b2))
				# printf "\nCalling differential peaks between ${sample_line_mark[i]} and ${sample_line_mark[j]}\n"
				# macs2 bdgdiff --t1 ChIP/peaks/${sample_line_mark[i]}_merged_treat_pileup.bdg --c1 ChIP/peaks/${sample_line_mark[i]}_merged_control_lambda.bdg --t2 ChIP/peaks/${sample_line_mark[j]}_merged_treat_pileup.bdg --c2 ChIP/peaks/${sample_line_mark[j]}_merged_control_lambda.bdg --d1 $a --d2 $b -g $g -l $l --outdir combined/peaks --o-prefix diff_${sample_line_mark[i]}_${sample_line_mark[j]}
				# j=$((j+1))
			# done
			# i=$((i+1))
		# done
	# else
		# printf "No differential ChIPseq peak analysis performed on $mark (not enough samples)\n"
	# fi
# done


#########################################################################################
####################################### PART6 ###########################################
################################### Making heatmaps  ####################################
#########################################################################################

#### To make heatmaps and profiles with deeptools
#### By default, it does both scale-regions and reference-point on start of bedfile provided
#### By default, it does heatmap on all the data, heatmap with 5 kmeans, and corresponding profiles
#### Probably need to edit many parameters depending on the purpose of the analysis

printf "\nDoing analysis for $analysisname with deeptools version:\n"
deeptools --version

#### To make heatmaps on the region file

#### Splitting the region file by strand
awk -v OFS="\t" '$6=="+"' $regionfile > combined/matrix/temp_regions_${regionname}_plus.bed
awk -v OFS="\t" '$6=="-"' $regionfile > combined/matrix/temp_regions_${regionname}_minus.bed

#### Reordering the samples by ChIPseq mark
for mark in ${uniq_chip_mark_list[@]}
do
	for sample in ${chip_sample_list[@]}
	do
		if [[ $sample =~ $mark ]]; then
			sorted_labels+=("$sample")
		fi
	done
	for bw in ${chip_bw_list[@]}
	do
		if [[ $bw =~ $mark ]]; then
			sorted_marks+=("$bw")
		fi
	done
done

## Computing the stranded matrix
for strand in plus minus
do
	case "$strand" in
		plus) 	bw_list="${sorted_marks[@]} ${rnaseq_bw_list_plus[@]}";;
		minus) 	bw_list="${sorted_marks[@]} ${rnaseq_bw_list_minus[@]}";;
	esac
	printf "\nComputing scale-regions $strand strand matrix for $analysisname\n"
	computeMatrix scale-regions --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${regionname}_${strand}.bed -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/regions_${analysisname}_${strand}.gz
	printf "\nComputing reference-point on TSS $strand strand matrix for $analysisname\n"
	computeMatrix reference-point --referencePoint "TSS" --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${regionname}_${strand}.bed -S ${bw_list} -bs 50 -b 2000 -a 8000 -p $threads -o combined/matrix/tss_${analysisname}_${strand}.gz
done
rm -f combined/matrix/temp_regions_${regionname}_*.bed

#### Merging stranded matrix, extracting scales and plotting heatmaps
for matrix in regions tss
do
	printf "\nMerging stranded matrices aligned by $matrix of $analysisname\n"
	computeMatrixOperations rbind -m combined/matrix/${matrix}_${analysisname}_plus.gz combined/matrix/${matrix}_${analysisname}_minus.gz -o combined/matrix/${matrix}_${analysisname}.gz
	printf "\nGetting scales (10th and 90th quantiles) for $matrix matrix of $analysisname\n"
	computeMatrixOperations dataRange -m combined/matrix/${matrix}_${analysisname}.gz > combined/matrix/values_${matrix}_${analysisname}.txt
	mins=()
	maxs=()
	if [ ${#rnaseq_bw_list_plus[@]} -gt 0 ]; then
		all_samples=("${uniq_chip_mark_list[*]}" "RNAseq")
		printf "\nincluding RNAseq to samples\n"
	else
		printf "\nOnly chipseq samples\n"
		all_samples=("${uniq_chip_mark_list[*]}")
	fi
	for mark in ${all_samples[@]}
	do
		mini=$(grep "$mark" combined/matrix/values_${matrix}_${analysisname}.txt | awk 'BEGIN {m=999999} {a=$5; if (a<m) m=a;} END {print m}')
		maxi=$(grep "$mark" combined/matrix/values_${matrix}_${analysisname}.txt | awk 'BEGIN {m=-999999} {a=$6; if (a>m) m=a;} END {print m}')
		num=$(grep "$mark" combined/matrix/values_${matrix}_${analysisname}.txt | wc -l)
		temp1=$(eval $(echo printf '"#$mini %0.s"' {1..$num}))
		temp2=$(eval $(echo printf '"#$maxi %0.s"' {1..$num}))
		temp3=${temp1[@]//#/}
		temp4=${temp2[@]//#/}
		mins+=("$temp3")
		maxs+=("$temp4")
	done
	printf "\nPlotting heatmap for $matrix matrix of $analysisname\n"
	plotHeatmap -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_${matrix}.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} ${rnaseq_sample_list[@]} --regionsLabel ${regionname} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --interpolationMethod 'bilinear'
	printf "\nPlotting heatmap for $matrix matrix of $analysisname in 5 clusters (kmeans)\n"
	plotHeatmap -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_${matrix}_k5.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} ${rnaseq_sample_list[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --interpolationMethod 'bilinear' --kmeans 5 --outFileSortedRegions combined/matrix/${analysisname}_${matrix}_regions_k5.txt
done

#### To make heatmaps and profiles with deeptools on the DEG if they were called


if [ ${#rnaseq_tissue_list[@]} -ge 2 ]; then
	#### To reorder bigwig files and sample names by ChIPseq mark
	sorted_marks=()
	sorted_labels=()
	for mark in ${uniq_chip_mark_list[@]}
	do
		for sample in ${chip_sample_list[@]}
		do
			if [[ $sample =~ $mark ]]; then
				sorted_labels+=("$sample")
			fi
		done
		for bw in ${chip_bw_list[@]}
		do
			if [[ $bw =~ $mark ]]; then
				sorted_marks+=("$bw")
			fi
		done
	done
	#### To plot all pairwise DEGs together
	regions_labels=()
	regions_files=()
	if [ -e combined/matrix/temp_regions_${analysisname}_all_DEGs.bed ]; then
		rm -f combined/matrix/temp_regions_${analysisname}_all_DEGs.bed
	fi
	for file in $(ls combined/DEG/DEG_${analysisname}*vs*.txt)
	do
		tmp3=${file##*${analysisname}_}
		filename=${tmp3%.txt}
		for DEG in UP DOWN
		do
			awk -v OFS="\t" -v d=$DEG '$11==d' $file > combined/matrix/temp_regions_${analysisname}_DEG_${filename}_${DEG}.bed
		done
		regions_files+=("combined/matrix/temp_regions_${analysisname}_DEG_${filename}_UP.bed" "combined/matrix/temp_regions_${analysisname}_DEG_${filename}_DOWN.bed")
		regions_labels+=("${filename}_UP" "${filename}_DOWN")
		awk -v OFS="\t" 'NR>1 {print $1,$2,$3}' $file >> combined/matrix/temp_regions_${analysisname}_all_DEGs.bed		
	done
	sort -k1,1n -k2,2n combined/matrix/temp_regions_${analysisname}_all_DEGs.bed -u > combined/matrix/temp_regions_${analysisname}_all_DEGs_unique.bed
	printf "\nComputing matrix for DEG for each sample pairs from $analysisname\n"
	computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${regions_files[@]} -S ${sorted_marks[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/${analysisname}_DEG.gz
	printf "\nComputing matrix for all DEGs from $analysisname\n"
	computeMatrix scale-regions --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${analysisname}_all_DEGs_unique.bed -S ${sorted_marks[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/${analysisname}_all_DEGs.gz
	printf "\nGetting scales (10th and 90th quantiles) for the DEG matrix of $analysisname\n"
	computeMatrixOperations dataRange -m combined/matrix/${analysisname}_DEG.gz > combined/matrix/values_${analysisname}_DEG.txt
	mins=()
	maxs=()
	for sample in ${sorted_labels[@]}
	do
		mini=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $5}')
		mins+=("$mini")
		maxi=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $6}')
		maxs+=("$maxi")
	done	
	printf "\nPlotting complete heatmap for DEG for each sample pairs from $analysisname\n"
	plotHeatmap -m combined/matrix/${analysisname}_DEG.gz -out combined/plots/${analysisname}_heatmap_DEG.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --interpolationMethod 'bilinear'
	plotHeatmap -m combined/matrix/${analysisname}_DEG.gz -out combined/plots/${analysisname}_heatmap_DEG_v2.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --interpolationMethod 'bilinear' --perGroup
	printf "\nPlotting complete profiles for DEG for each sample pairs from $analysisname\n"
	plotProfile -m combined/matrix/${analysisname}_DEG.gz -out combined/plots/${analysisname}_profile_DEG.pdf --plotType 'lines' --averageType 'median' --samplesLabel ${sorted_labels[@]} --regionsLabel ${regions_labels[@]}
	plotProfile -m combined/matrix/${analysisname}_DEG.gz -out combined/plots/${analysisname}_profile_DEG_v2.pdf --plotType 'lines' --averageType 'median' --samplesLabel ${sorted_labels[@]} --regionsLabel ${regions_labels[@]} --perGroup
	printf "\nPlotting complete heatmap for DEG for each sample pairs from $analysisname\n"
	plotHeatmap -m combined/matrix/${analysisname}_all_DEGs.gz -out combined/plots/${analysisname}_heatmap_all_DEGs_k5.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} --colorMap 'seismic' --interpolationMethod 'bilinear' --kmeans 5
	printf "\nPlotting complete profiles for DEG for each sample pairs from $analysisname\n"
	plotProfile -m combined/matrix/${analysisname}_all_DEGs.gz -out combined/plots/${analysisname}_profile_all_DEGs_k5.pdf --plotType 'lines' --averageType 'median' --samplesLabel ${sorted_labels[@]} --kmeans 5
	
	# #### To plot tissue-specific DEGs
	
	# for file in $(ls combined/DEG/only_${tissue}_DEG_${analysisname}.bed)
	# do
		# tmp4=${file%_DEG_${analysisname}.bed}
		# filename=${tmp4##*/}
		# tissue=${tmp4##*/only_}
		# printf "\nComputing matrix for $tissue only ($filename) DEG from $analysisname\n"
		# computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${file} -S ${sorted_marks[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/${analysisname}_${filename}_DEG.gz
		# printf "\nPlotting heatmap for $tissue only ($filename) DEG from $analysisname\n"
		# plotHeatmap -m combined/matrix/${analysisname}_${filename}_DEG.gz -out combined/plots/${analysisname}_heatmap_${filename}_DEG.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} --regionsLabel ${filename} --colorMap 'seismic' --interpolationMethod 'bilinear'
	# done
	
fi
rm -f combined/matrix/temp_regions_${analysisname}_*.bed

printf "\nCombined analysis script finished successfully for $analysisname\n"
touch combined/chkpts/analysis_${analysisname}
