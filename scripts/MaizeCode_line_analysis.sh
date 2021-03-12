#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=4G
#$ -l tmp_free=100G
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
##### Requirements: bedtools, deeptools, macs2, R (+R packages: ggplot2,UpSetR,limma,edgeR,dplyr,tidyr,stringr,gplots)
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS
export mc_dir="${PWD}/scripts/"
printf "\nRunning MaizeCode scripts in working directory ${PWD}\n"

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
rnaseq_name_list=()
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
			rnaseq_name_list+=("${line}_${tissue}")
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
	awk -v OFS="\t" 'BEGIN {printf "PeakID\tDistance\tGroup\n"} {if ($5<-2000) {d="Distal_upstream"; a=-$5}; else if ($5<0) {d="Promoter"; a=-$5}; else if ($5==0) {d="Gene_body"; a=$5}; else if ($5>2000) {d="Distal_downstream"; a=$5}; else {d="Terminator"; a=$5}; print $4,a,d}' combined/peaks/peaks_${analysisname}.bed > combined/peaks/temp_col_${analysisname}_AAA.txt
	paste combined/peaks/temp_col_${analysisname}_*.txt | uniq > combined/peaks/matrix_upset_${analysisname}.txt
	rm -f combined/peaks/temp_col_${analysisname}_*.txt
	#### To make an Upset plot highlighting peaks in gene bodies
	printf "\nCreating Upset plot for $samplename with R version:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_Upset.r combined/peaks/matrix_upset_${analysisname}.txt ${analysisname}
	Rscript --vanilla ${mc_dir}/MaizeCode_R_Upset_v2.r combined/peaks/matrix_upset_${analysisname}.txt ${analysisname}
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

if [ ${#rnaseq_sample_list[@]} -ge 2 ]; then
	#### To make a count table for all RNAseq samples in samplefile
	printf "\nPreparing count table for RNAseq samples in $analysisname\n"
	printf "Replicate\tSample\tColor\n" > combined/DEG/samples_${analysisname}.txt
	if [ $(grep "gene:" RNA/mapped/map_${rnaseq_sample_list[0]}_Rep1_ReadsPerGene.out.tab | wc -l) -gt 0 ]; then
		i=0
		for sample in ${rnaseq_sample_list[@]}
		do
			namei=${rnaseq_name_list[i]}
			numreps=$(ls -1f RNA/mapped/map_${sample}_Rep*_ReadsPerGene.out.tab | wc -l)
			n=$((i+1))
			for ((j=1;j<=numreps;j++))
			do
				printf "${namei}_Rep${j}\t${namei}\t${n}\n" >> combined/DEG/samples_${analysisname}.txt
				grep "gene:" RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab | sed 's/gene://' | awk -v OFS="\t" -v t=$namei -v j=$j 'BEGIN {print t"_Rep"j} {print $2}' > combined/DEG/col_AZ_${i}_${analysisname}_${sample}_Rep${j}.txt
				if [ $i -eq 0 ] && [ $j -eq 1 ]; then
					grep "gene:" RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab | sed 's/gene://' | awk -v OFS="\t" 'BEGIN {print "gene_ID"} {print $1}' > combined/DEG/col_AA_0_${analysisname}.txt
				fi
			done
			i=$((i+1))
		done
	else
		i=0
		for sample in ${rnaseq_sample_list[@]}
		do
			namei=${rnaseq_name_list[i]}
			numreps=$(ls -1f RNA/mapped/map_${sample}_Rep*_ReadsPerGene.out.tab | wc -l)
			n=$((i+1))
			for ((j=1;j<=numreps;j++))
			do
				printf "${namei}_Rep${j}\t${namei}\t${n}\n" >> combined/DEG/samples_${analysisname}.txt
				awk -v OFS="\t" -v t=$namei -v j=$j 'BEGIN {print t"_Rep"j}  $1 !~ /^N_/ {print $2}' RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab > combined/DEG/col_AZ_${i}_${analysisname}_${sample}_Rep${j}.txt
				if [ $i -eq 0 ] && [ $j -eq 1 ]; then
					awk -v OFS="\t" 'BEGIN {print "gene_ID"}  $1 !~ /^N_/ {print $1}' RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab > combined/DEG/col_AA_0_${analysisname}.txt
				fi
			done
			i=$((i+1))
		done
	fi
	paste combined/DEG/col_A*_${analysisname}* > combined/DEG/counts_${analysisname}.txt
	rm -f combined/DEG/col_A*_${analysisname}*
	#### To run the DEG analysis on R
	printf "\nLaunching DEG analysis with R version:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_DEG.r combined/DEG/counts_${analysisname}.txt combined/DEG/samples_${analysisname}.txt ${analysisname} $regionfile
	#### To extract DEG only called in one tissue
	if [ ${#rnaseq_name_list[@]} -ge 3 ]; then
		for namei in ${rnaseq_name_list[@]}
		do
			if ls combined/DEG/DEG_${analysisname}_${namei}_vs*.txt 1> /dev/null 2>&1; then
				for file in $(ls combined/DEG/DEG_${analysisname}_${namei}_vs*.txt)
				do
					awk -v OFS="\t" '$11 == "DOWN" {print $1,$2,$3,$4,".",$6}' ${file} > ${file}.DOWN.temp.bed
					awk -v OFS="\t" '$11 == "UP" {print $1,$2,$3,$4,".",$6}' ${file} > ${file}.UP.temp.bed
				done
			fi
			if ls combined/DEG/DEG_${analysisname}_*_vs_${namei}.txt 1> /dev/null 2>&1; then
				for file in $(ls combined/DEG/DEG_${analysisname}_*vs_${namei}.txt)
				do
					awk -v OFS="\t" '$11 == "DOWN" {print $1,$2,$3,$4,".",$6}' ${file} >> ${file}.UP.temp.bed
					awk -v OFS="\t" '$11 == "UP" {print $1,$2,$3,$4,".",$6}' ${file} >> ${file}.DOWN.temp.bed
				done
			fi
			printf "\nGetting DEGs specific in $namei\n"
			filearraydown=( $(ls combined/DEG/DEG_${analysisname}_*vs*.txt.DOWN.temp.bed | grep "$namei") )
			filearrayup=( $(ls combined/DEG/DEG_${analysisname}_*vs*.txt.UP.temp.bed | grep "$namei") )
			i=0
			max=$((${#filearraydown[@]}-1))
			cat ${filearraydown[i]} > combined/DEG/temp_tissue_spec_DEG_${analysisname}_DOWN_${i}.txt
			cat ${filearrayup[i]} > combined/DEG/temp_tissue_spec_DEG_${analysisname}_UP_${i}.txt
			while [ $i -lt $max ]
			do
				j=$((i+1))
				bedtools intersect -wa -a combined/DEG/temp_tissue_spec_DEG_${analysisname}_DOWN_${i}.txt -b ${filearraydown[j]} > combined/DEG/temp_tissue_spec_DEG_${analysisname}_DOWN_${j}.txt
				bedtools intersect -wa -a combined/DEG/temp_tissue_spec_DEG_${analysisname}_UP_${i}.txt -b ${filearrayup[j]} > combined/DEG/temp_tissue_spec_DEG_${analysisname}_UP_${j}.txt
				i=$((i+1))
			done
			cat combined/DEG/temp_tissue_spec_DEG_${analysisname}_DOWN_${max}.txt > combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed
			cat combined/DEG/temp_tissue_spec_DEG_${analysisname}_UP_${max}.txt > combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed
			rm -f combined/DEG/DEG_${analysisname}_*.temp.bed
			rm -f combined/DEG/temp_tissue_spec_DEG_${analysisname}*
		done
	fi
else 
	printf "\nNo differential gene expression analysis performed (not enough samples)\n"
fi


############################################################################################
########################################## PART5 ###########################################
##################### Differential peaks analysis between tissues ##########################
############################################################################################

# #### Takes a long time and not used yet so commented until further improvements ####

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

uniq_chip_mark_list=($(printf "%s\n" "${chip_mark_list[@]}" | sort -u))

#### To make heatmaps on the region file

#### Splitting the region file by strand
#awk -v OFS="\t" '$6=="+"' $regionfile > combined/matrix/temp_regions_${regionname}_plus.bed
#awk -v OFS="\t" '$6=="-"' $regionfile > combined/matrix/temp_regions_${regionname}_minus.bed

#### Reordering the samples by ChIPseq mark
#for mark in ${uniq_chip_mark_list[@]}
#do
#	for sample in ${chip_sample_list[@]}
#	do
#		if [[ $sample =~ $mark ]]; then
#			sorted_labels+=("$sample")
#		fi
#	done
#	for bw in ${chip_bw_list[@]}
#	do
#		if [[ $bw =~ $mark ]]; then
#			sorted_marks+=("$bw")
#		fi
#	done
#done

#### Computing the stranded matrix
#for strand in plus minus
#do
#	case "$strand" in
#		plus) 	bw_list="${sorted_marks[@]} ${rnaseq_bw_list_plus[@]}";;
#		minus) 	bw_list="${sorted_marks[@]} ${rnaseq_bw_list_minus[@]}";;
#	esac
#	printf "\nComputing scale-regions $strand strand matrix for $analysisname\n"
#	computeMatrix scale-regions --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${regionname}_${strand}.bed -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/regions_${analysisname}_${strand}.gz
#	printf "\nComputing reference-point on TSS $strand strand matrix for $analysisname\n"
#	computeMatrix reference-point --referencePoint "TSS" --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${regionname}_${strand}.bed -S ${bw_list} -bs 50 -b 2000 -a 8000 -p $threads -o combined/matrix/tss_${analysisname}_${strand}.gz
#done
#rm -f combined/matrix/temp_regions_${regionname}_*.bed

#### Merging stranded matrix, extracting scales and plotting heatmaps
#for matrix in regions tss
#do
#	printf "\nMerging stranded matrices aligned by $matrix of $analysisname\n"
#	computeMatrixOperations rbind -m combined/matrix/${matrix}_${analysisname}_plus.gz combined/matrix/${matrix}_${analysisname}_minus.gz -o combined/matrix/${matrix}_${analysisname}.gz
#	printf "\nGetting scales (10th and 90th quantiles) for $matrix matrix of $analysisname\n"
#	computeMatrixOperations dataRange -m combined/matrix/${matrix}_${analysisname}.gz > combined/matrix/values_${matrix}_${analysisname}.txt
#	mins=()
#	maxs=()
#	if [ ${#rnaseq_bw_list_plus[@]} -gt 0 ]; then
#		all_samples=("${uniq_chip_mark_list[*]}" "RNAseq")
#		printf "\nincluding RNAseq to samples\n"
#	else
#		printf "\nOnly chipseq samples\n"
#		all_samples=("${uniq_chip_mark_list[*]}")
#	fi
#	for mark in ${all_samples[@]}
#	do
#		mini=$(grep "$mark" combined/matrix/values_${matrix}_${analysisname}.txt | awk 'BEGIN {m=999999} {a=$5; if (a<m) m=a;} END {print m}')
#		maxi=$(grep "$mark" combined/matrix/values_${matrix}_${analysisname}.txt | awk 'BEGIN {m=-999999} {a=$6; if (a>m) m=a;} END {print m}')
#		num=$(grep "$mark" combined/matrix/values_${matrix}_${analysisname}.txt | wc -l)
#		temp1=$(eval $(echo printf '"#$mini %0.s"' {1..$num}))
#		temp2=$(eval $(echo printf '"#$maxi %0.s"' {1..$num}))
#		temp3=${temp1[@]//#/}
#		temp4=${temp2[@]//#/}
#		mins+=("$temp3")
#		maxs+=("$temp4")
#	done
#	mins2=()
#	maxs2=()
#	for sample in ${sorted_labels[@]} ${rnaseq_sample_list[@]}
#	do
#		mini=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $5}')
#		mins2+=("$mini")
#		maxi=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $6}')
#		maxs2+=("$maxi")
#	done
#	printf "\nPlotting heatmap for $matrix matrix of $analysisname scaling by mark\n"
#	plotHeatmap -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_${matrix}.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} ${rnaseq_sample_list[@]} --regionsLabel ${regionname} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --interpolationMethod 'bilinear'
#	printf "\nPlotting heatmap for $matrix matrix of $analysisname scaling by sample\n"
#	plotHeatmap -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_${matrix}_v2.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} ${rnaseq_sample_list[@]} --regionsLabel ${regionname} --colorMap 'seismic' --zMin ${mins2[@]} --zMax ${maxs2[@]} --interpolationMethod 'bilinear'
	# #### Not very useful right now ####
	# printf "\nPlotting heatmap for $matrix matrix of $analysisname in 5 clusters (kmeans)\n"
	# plotHeatmap -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_${matrix}_k5.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} ${rnaseq_sample_list[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --interpolationMethod 'bilinear' --kmeans 5 --outFileSortedRegions combined/matrix/${analysisname}_${matrix}_regions_k5.txt
#done

#rm -f combined/matrix/*${analysisname}*.gz

#### To make heatmaps and profiles with deeptools on the DEG if they were called

#if [ ${#rnaseq_name_list[@]} -ge 2 ]; then
	#### To reorder bigwig files and sample names by ChIPseq mark
#	sorted_marks=()
#	sorted_labels=()
#	for mark in ${uniq_chip_mark_list[@]}
#	do
#		for sample in ${chip_sample_list[@]}
#		do
#			if [[ $sample =~ $mark ]]; then
#				sorted_labels+=("$sample")
#			fi
#		done
#		for bw in ${chip_bw_list[@]}
#		do
#			if [[ $bw =~ $mark ]]; then
#				sorted_marks+=("$bw")
#			fi
#		done
#	done
#	if [ ${#sorted_marks[@]} -ge 1 ]; then
#		#### To plot all pairwise DEGs together
#		regions_labels=()
#		regions_files=()
#		if [ -e combined/matrix/temp_regions_${analysisname}_all_DEGs.bed ]; then
#			rm -f combined/matrix/temp_regions_${analysisname}_all_DEGs.bed
#		fi
#		for file in $(ls combined/DEG/DEG_${analysisname}*vs*.txt)
#		do
#			tmp3=${file##*${analysisname}_}
#			filename=${tmp3%.txt}
#			for DEG in UP DOWN
#			do
#				awk -v OFS="\t" -v d=$DEG '$11==d' $file > combined/matrix/temp_regions_${analysisname}_DEG_${filename}_${DEG}.bed
#			done
#			regions_files+=("combined/matrix/temp_regions_${analysisname}_DEG_${filename}_UP.bed" "combined/matrix/temp_regions_${analysisname}_DEG_${filename}_DOWN.bed")
#			regions_labels+=("${filename}_UP" "${filename}_DOWN")
#			awk -v OFS="\t" 'NR>1 {print $1,$2,$3}' $file >> combined/matrix/temp_regions_${analysisname}_all_DEGs.bed		
#		done
#		sort -k1,1n -k2,2n combined/matrix/temp_regions_${analysisname}_all_DEGs.bed -u > combined/matrix/temp_regions_${analysisname}_all_DEGs_unique.bed
#		printf "\nComputing matrix for DEG for each sample pairs from $analysisname\n"
#		computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${regions_files[@]} -S ${sorted_marks[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/${analysisname}_DEG.gz
#		printf "\nComputing matrix for all DEGs from $analysisname\n"
#		computeMatrix scale-regions --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${analysisname}_all_DEGs_unique.bed -S ${sorted_marks[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/${analysisname}_all_DEGs.gz
#		printf "\nGetting scales (10th and 90th quantiles) for the DEG matrix of $analysisname\n"
#		computeMatrixOperations dataRange -m combined/matrix/${analysisname}_DEG.gz > combined/matrix/values_${analysisname}_DEG.txt
#		mins=()
#		maxs=()
#		for sample in ${sorted_labels[@]}
#		do
#			mini=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $5}')
#			mins+=("$mini")
#			maxi=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $6}')
#			maxs+=("$maxi")
#		done
		# #### Not very useful right now, too messy ####
		# k=$(( ${#regions_files[@]} / 2 ))
		# printf "\nPlotting complete heatmap for DEG for each sample pairs from $analysisname\n"
		# plotHeatmap -m combined/matrix/${analysisname}_DEG.gz -out combined/plots/${analysisname}_heatmap_DEG.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --interpolationMethod 'bilinear'
		# printf "\nPlotting complete heatmap for all DEGs from $analysisname\n"
		# plotHeatmap -m combined/matrix/${analysisname}_all_DEGs.gz -out combined/plots/${analysisname}_heatmap_all_DEGs_${k}.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} --colorMap 'seismic' --interpolationMethod 'bilinear' --kmeans ${k}
	
#		for mark in ${uniq_chip_mark_list[@]}
#		do
#			selected_samples=()
#			selected_labels=()
#			for sample in ${sorted_labels[@]}
#			do
#				if [[ $sample =~ $mark ]]; then
#					selected_samples+=("${sample}_merged")
#					selected_labels+=("${sample}")
#				fi
#			done		
#			computeMatrixOperations subset -m combined/matrix/${analysisname}_DEG.gz -o combined/matrix/${analysisname}_DEG_${mark}.gz --samples ${selected_samples[@]}
#			printf "\nPlotting ${mark} profiles for DEG for each sample pairs from $analysisname\n"
#			plotProfile -m combined/matrix/${analysisname}_DEG_${mark}.gz -out combined/plots/${analysisname}_profile_DEG_${mark}.pdf --plotType 'lines' --averageType 'median' --samplesLabel ${selected_labels[@]} --regionsLabel ${regions_labels[@]} --perGroup --numPlotsPerRow 2
#		done
#		#### To plot tissue-specific DEGs
#		if [ ${#rnaseq_name_list[@]} -ge 3 ]; then
#			for namei in ${rnaseq_name_list[@]}
#			do
#				filenames="combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed"
#				printf "\nComputing matrix for $namei specific DEG from $analysisname\n"
#				computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${filenames} -S ${sorted_marks[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/${analysisname}_only_${namei}_DEG.gz
#				printf "\nGetting scales (10th and 90th quantiles) for the $tissue specific DEG matrix of $analysisname\n"
#				computeMatrixOperations dataRange -m combined/matrix/${analysisname}_only_${namei}_DEG.gz > combined/matrix/values_${analysisname}_only_${namei}_DEG.gz
#				mins=()
#				maxs=()
#				for sample in ${sorted_labels[@]}
#				do
#					mini=$(grep $sample combined/matrix/values_${analysisname}_only_${namei}_DEG.gz | awk '{print $5}')
#					mins+=("$mini")
#					maxi=$(grep $sample combined/matrix/values_${analysisname}_only_${namei}_DEG.gz | awk '{print $6}')
#					maxs+=("$maxi")
#				done	
#				printf "\nPlotting heatmap for $namei specific DEG from $analysisname\n"
#				plotHeatmap -m combined/matrix/${analysisname}_only_${namei}_DEG.gz -out combined/plots/${analysisname}_heatmap_only_${namei}_DEG.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} --regionsLabel "${namei}_UP" "${namei}_DOWN" --zMin ${mins[@]} --zMax ${maxs[@]} --colorMap 'seismic' --interpolationMethod 'bilinear'
#			done
#		fi
#	fi
#fi
#rm -f combined/matrix/temp_regions_${analysisname}*.bed
#rm -f combined/matrix/*${analysisname}*.gz
#rm -f combined/matrix/values*${analysisname}*

#### To make heatmap and profile with deeptools for each tissue based on grouped expression levels (if both RNA and ChIP samples are present in a tissue)
#### Only works when chromosomes starts with 1-10, not Chr. Need to be updated.

if [[ $ref == "B73_v4" ]]; then
	for tissue in ${rnaseq_tissue_list[*]}
	do
		if [[ " ${chip_tissue_list[@]} " =~ " ${tissue} " ]]; then
			tissue_chip_bw_list=()
			for bw in ${chip_bw_list[*]}
			do
				if [[ $bw =~ $tissue ]]; then
					tissue_chip_bw_list+=("$bw")
				fi
			done
			tissue_labels_chip=()
			for sample in ${chip_sample_list[*]}
			do
				if [[ $sample =~ $tissue ]]; then
					tissue_labels_chip+=("$sample")
				fi
			done
			tissue_rnaseq_bw_list_plus=()
			for bw in ${rnaseq_bw_list_plus[*]}
			do					
				if [[ $bw =~ $tissue ]]; then
					tissue_rnaseq_bw_list_plus+=("$bw")
				fi
			done
			tissue_rnaseq_bw_list_minus=()
			for bw in ${rnaseq_bw_list_minus[*]}
			do					
				if [[ $bw =~ $tissue ]]; then
					tissue_rnaseq_bw_list_minus+=("$bw")
				fi
			done
			tissue_labels_rna=()
			for sample in ${rnaseq_sample_list[*]}
			do
				if [[ $sample =~ $tissue ]]; then
					tissue_labels_rna+=("$sample")
				fi
			done

			printf "Clustering genes by expression levels for ${tissue}\n"
			cols=($(awk -v ORS=" " -v t=$tissue 'NR==1 {for(i=1;i<=NF;i++) if ($i~t) print i}' combined/DEG/counts_${analysisname}.txt))
			reps=${#cols[@]}
			awk -v d="$cols" -v t=$reps 'BEGIN {split(d, a, " ")} NR > 1 {b=0; for (i in a) b+=$(a[i]); c=b/t; print $1,c}' combined/DEG/counts_${analysisname}.txt > combined/DEG/temp_counts_${analysisname}_${tissue}.txt
			if [ -s combined/DEG/temp_expression_${analysisname}_${tissue}.bed ]; then
				rm -f combined/DEG/temp_expression_${analysisname}_${tissue}.bed
			fi
			while read ID exp
			do
				grep "$ID" $regionfile | awk -v OFS="\t" -v c=$exp '$1 ~ /^[0-9]/ {l=$3-$2; $5=1000*c/l; print $0}' >> combined/DEG/temp_expression_${analysisname}_${tissue}.bed
			done < combined/DEG/temp_counts_${analysisname}_${tissue}.txt
			sort -k5,5gr combined/DEG/temp_expression_${analysisname}_${tissue}.bed > combined/DEG/sorted_expression_${analysisname}_${tissue}.bed
		
			awk -v OFS="\t" -v a=${analysisname} -v b=${tissue} '{if ($5==0) printf $0"\n" > "combined/DEG/sorted_"a"_"b"_exp0.bed"; else printf $0"\n" > "combined/DEG/sorted_"a"_"b"_expA.bed"}' combined/DEG/sorted_expression_${analysisname}_${tissue}.bed
			tot=$(wc -l combined/DEG/sorted_${analysisname}_${tissue}_expA.bed | awk '{print $1}')
			bin=$((tot/5))
			min=0
			max=$bin
			for (( i = 1; i <= 5; i++ ))
			do
				awk -v n=$min -v m=$max 'NR>=n && NR <=m' combined/DEG/sorted_${analysisname}_${tissue}_expA.bed > combined/DEG/sorted_${analysisname}_${tissue}_exp${i}.bed
				min=$((min+bin))
				max=$((max+bin))
			done

			sorted_regions=()
			regions_labels=()
			for i in 1 2 3 4 5 0
			do
				case "$i" in
					1) name="Top20%";;
					2) name="20-40%";;
					3) name="40-60%";;
					4) name="60-80%";;
					5) name="Bottom20%";;
					0) name="Not_expressed";;
				esac
				n=$(wc -l combined/DEG/sorted_${analysisname}_${tissue}_exp${i}.bed | awk '{print $1}')
				regions_labels+=("$name($n)")	
				sorted_regions+=("combined/DEG/sorted_${analysisname}_${tissue}_exp${i}.bed")
				\cp -Tf combined/DEG/sorted_${analysisname}_${tissue}_exp${i}.bed combined/DEG/sorted_${analysisname}_${tissue}_exp${i}.txt
			done	
			for strand in plus minus
			do
				case "$strand" in
					plus) 	bw_list="${tissue_chip_bw_list[*]} ${tissue_rnaseq_bw_list_plus[*]}"
						sign="+";;
					minus) 	bw_list="${tissue_chip_bw_list[*]} ${tissue_rnaseq_bw_list_minus[*]}"
						sign="-";;
				esac
				regions=()
				for i in 1 2 3 4 5 0
				do
					awk -v OFS="\t" -v s=$sign '$6==s' combined/DEG/sorted_${analysisname}_${tissue}_exp${i}.txt > combined/DEG/sorted_${analysisname}_${tissue}_exp${i}.bed
				done
				printf "\nComputing scale-regions $strand strand matrix for ${tissue} in ${analysisname}\n"
				computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/regions_${analysisname}_${strand}.gz
				printf "\nComputing reference-point on TSS $strand strand matrix for ${tissue} in $analysisname\n"
				computeMatrix reference-point --referencePoint "TSS" --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 8000 -p $threads -o combined/matrix/tss_${analysisname}_${strand}.gz
			done
		
			for i in 1 2 3 4 5 0
			do
				\cp -Tf combined/DEG/sorted_${analysisname}_${tissue}_exp${i}.txt combined/DEG/sorted_${analysisname}_${tissue}_exp${i}.bed
			done

			### Merging stranded matrix, extracting scales and plotting heatmaps
			for matrix in regions tss
			do
				tissue_labels="${tissue_labels_chip[*]} ${tissue_labels_rna[*]}"
				printf "\nMerging stranded matrices aligned by $matrix for ${tissue} in $analysisname\n"
				computeMatrixOperations rbind -m combined/matrix/${matrix}_${analysisname}_plus.gz combined/matrix/${matrix}_${analysisname}_minus.gz -o combined/matrix/${matrix}_${analysisname}.gz
				printf "\nGetting scales (10th and 90th quantiles) for $matrix matrix for ${tissue} in $analysisname\n"
				computeMatrixOperations dataRange -m combined/matrix/${matrix}_${analysisname}.gz > combined/matrix/values_${matrix}_${analysisname}.txt
				mins=()
				maxs=()
				for sample in ${tissue_labels[@]}
				do
					mini=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $5}')
					mins+=("$mini")
					maxi=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $6}')
					maxs+=("$maxi")
				done
				computeMatrixOperations sort -m combined/matrix/${matrix}_${analysisname}.gz -R ${sorted_regions[@]} -o combined/matrix/final_${matrix}_${analysisname}.gz
				plotProfile -m combined/matrix/final_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_profile_${matrix}.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --outFileNameData combined/matrix/values_${matrix}_${analysisname}.txt
				ymins=()
				ymax=()
				for sample in ${tissue_labels[@]}
				do
				 	ymini=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m*1.1}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {print m}')
					ymins+=("$ymini")
					ymaxi=$(grep $sample combined/matrix/values_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m*1.1}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m}')
					ymaxs+=("$ymaxi")		
				done
				printf "yMins = ${ymins[*]}\nyMaxs = ${ymaxs[*]}\n"
				printf "\nPlotting heatmap for $matrix matrix for ${tissue} in $analysisname scaling by sample\n"
				plotHeatmap -m combined/matrix/final_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_heatmap_${matrix}.pdf --sortRegions keep --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
				printf "\nPlotting profile for $matrix matrix for ${tissue} in $analysisname scaling by sample\n"
				plotProfile -m combined/matrix/final_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_profile_${matrix}.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --yMin ${ymins[@]} --yMax ${ymaxs[@]}
			done
		fi
	done
fi

#### To make heatmaps and profiles with deeptools on the ACRs called in Ricci et al. 2019 paper (if B73_v4 is the reference used)
#### Not yet used so commented

# if [[ $ref =~ "B73_v4" ]]; then
#	if [ ! -s combined/matrix/leaf_ACRs.bed ]; then
#		wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3398nnn/GSM3398046/suppl/GSM3398046_ATAC_B73_leaf.filtered_ACR.bed.gz
#		pigz -d GSM3398046_ATAC_B73_leaf.filtered_ACR.bed.gz && mv GSM3398046_ATAC_B73_leaf.filtered_ACR.bed combined/matrix/leaf_ACRs.bed
#	fi
#	if [ ! -s combined/matrix/ears_ACRs.bed ]; then
#		wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3398nnn/GSM3398046/suppl/GSM3398047_ATAC_B73_ear.filtered_ACR.bed.gz
#		pigz -d GSM3398047_ATAC_B73_ear.filtered_ACR.bed.gz && mv GSM3398047_ATAC_B73_ear.filtered_ACR.bed combined/matrix/ears_ACRs.bed
#	fi
# fi	


printf "\nCombined analysis script finished successfully for $analysisname\n"
touch combined/chkpts/analysis_${analysisname}
