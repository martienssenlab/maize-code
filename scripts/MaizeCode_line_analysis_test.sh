#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=4G
#$ -l tmp_free=8G
#$ -o lineanalysis.log
#$ -j y
#$ -N lineanalysis

usage="
##### Script for Maize code inbred line data analysis, used by script MaizeCode_analyis.sh if -s was not set and region file exists
#####
##### sh MaiCode_line_analysis.sh -f samplefile -r regionfile [-t]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Sample, PE or SE, Reference genome directory
##### 	-r: bedfile containing the regions that are to be ploted over
#####	-t: If set, partial analysis will be performed (no heatmap with deeptools)	
##### 	-h: help, returns usage
##### 
##### After preparing all the data in the samplefile (PART1) it will, depending on the types of data present:
##### PART2: perform differential gene expression analysis among all RNAseq samples, potentially with gene ontology analysis (for B73_v4 for now)
##### PART3: produce an Upset plot of the intersection between all ChIP samples, highlighting peaks in the input regions
##### PART4: produce an Upset plot of the intersection between all TF samples, highlighting peaks in the H3K27ac peaks (if possible)
##### PART8: plot different heatmaps and metaplots of all ChIP and RNA samples in the samplefile on all input regions
##### PART9: plot heatmaps and metaplots on DEGs
##### PART10: plot heatmaps and metaplots on genes split by expression levels
##### PART11: plot heatmaps and metaplots on distal peaks split by H3K27ac ChIPseq quality
##### Under development:
##### PART5: compare TF peaks with DEGs
##### PART6: perform differential peak/TSS calling between all pairs of RAMPAGE samples
##### PART7: perform differential peak calling between all pairs of ChIP samples for each mark
##### PART11: perform differential sRNA expression analysis
#####
##### Requirements: bedtools, deeptools, macs2, R (+R packages: ggplot2,UpSetR,limma,edgeR,dplyr,tidyr,stringr,gplots)
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=${NSLOTS}
export mc_dir="${PWD}/scripts/"
printf "\nRunning MaizeCode scripts in working directory ${PWD}\n"

if [ $# -eq 0 ]; then
	printf "${usage}\n"
	exit 1
fi

while getopts ":f:r:th" opt; do
	case ${opt} in
		h) 	printf "${usage}\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		r)	export regionfile=${OPTARG};;
		t)	printf "\nPartial analysis to be performed (no heatmaps)\n" 
			export total="No";;
		*)	printf "\nUnknown argument!\n${usage}\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! ${samplefile} ]; then
	printf "Missing samplefile!\n"
	printf "${usage}\n"
	exit 1
fi

if [ ! ${regionfile} ]; then
	printf "Missing regionfile!\n"
	printf "${usage}\n"
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

printf "\nStarting analysis: ${analysisname}\n"

#### To append the lists (used for labels and other stuff)

ref_dir_list=()
line_list=()
chip_sample_list=()
chip_tissue_list=()
chip_mark_list=()
chip_bw_list=()
tf_sample_list=()
tf_tissue_list=()
tf_bw_list=()
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
rampage_bw_list_plus=()
rampage_bw_list_minus=()
shrna_bw_list_plus=()
shrna_bw_list_minus=()
while read data line tissue sample paired ref_dir
do
	case "${data}" in
		ChIP) 	datatype="ChIP"
			name=${line}_${tissue}_${sample};;
		RNAseq) datatype="RNA"
			name=${line}_${tissue}_${sample};;
		RAMPAGE) datatype="RNA"
			name=${line}_${tissue}_${sample};;
		shRNA) datatype="shRNA"
			name=${line}_${tissue}_${sample};;
		TF_*) datatype="TF"
			tmpname=${data##TF_}
			name=${line}_${tmpname};;
	esac

	ref=${ref_dir##*/}
	if [[ ! "${ref_list[@]}" =~ "${ref}" ]]; then
		ref_list+=("${ref}")
		ref_dir_list+=("${ref_dir}")
	fi
	if [[ ! "${line_list[@]}" =~ "${line}" ]]; then
		line_list+=("${line}")
	fi
	if [[ "${datatype}" == "ChIP" ]]; then
		chip_bw_list+=("${datatype}/tracks/${name}_merged.bw")
		chip_sample_list+=("${name}")
		chip_tissue_list+=("${tissue}")
		chip_mark_list+=("${sample}")
	elif [[ "${datatype}" == "RNA" ]]; then
		rna_type_list+=("${sample}")
		if [[ "${sample}" == "RNAseq" ]]; then
			rnaseq_bw_list_plus+=("${datatype}/tracks/${name}_merged_plus.bw")
			rnaseq_bw_list_minus+=("${datatype}/tracks/${name}_merged_minus.bw")
			rnaseq_sample_list+=("${name}")
			rnaseq_tissue_list+=("${tissue}")
			rnaseq_name_list+=("${line}_${tissue}")
		elif [[ "${sample}" == "RAMPAGE" ]]; then
			rampage_bw_list_plus+=("${datatype}/tracks/${name}_merged_plus.bw")
			rampage_bw_list_minus+=("${datatype}/tracks/${name}_merged_minus.bw")
			rampage_sample_list+=("${name}")
			rampage_tissue_list+=("${tissue}")
		else 
			printf "\nType of RNA sample unknown\n"
		fi
	elif [[ "${datatype}" == "shRNA" ]]; then
		shrna_sample_list+=("${name}")
		shrna_tissue_list+=("${tissue}")
		shrna_bw_list_plus+=("${datatype}/tracks/${name}_merged_plus.bw")
		shrna_bw_list_minus+=("${datatype}/tracks/${name}_merged_minus.bw")
	elif [[ "${datatype}" == "TF" ]]; then
		tf_sample_list+=("${name}")
		tf_tissue_list+=("${tissue}")
		tf_bw_list+=("${datatype}/tracks/${name}_merged.bw")
	else
		printf "\nType of data unknown for ${name}\nSample not processed!\n"
	fi
done < ${samplefile}

if [ ! ${#ref_list[@]} -eq 1 ] || [ ! ${#line_list[@]} -eq 1 ]; then
	printf "\nThere are multiple references in the samplefile! This analysis cannot be performed!\n"
	exit 1
else
	export ref=${ref_list[0]}
	export ref_dir=${ref_dir_list[0]}
	export line=${line_list[0]}
fi

#################

uniq_chip_mark_list=($(printf "%s\n" "${chip_mark_list[@]}" | sort -u))

for mark in ${uniq_chip_mark_list[@]}
do
	for sample in ${chip_sample_list[@]}
	do
		if [[ ${sample} =~ ${mark} ]]; then
			sorted_labels+=("${sample}")
		fi
	done
	for bw in ${chip_bw_list[@]}
	do
		if [[ ${bw} =~ ${mark} ]]; then
			sorted_marks+=("${bw}")
		fi
	done
done

#########################################################################################
####################################### PART9 ###########################################
############################## Making heatmaps on DEGs  #################################
#########################################################################################

#### To make heatmaps and profiles with deeptools on the DEG if they were called

if [ ${#rnaseq_name_list[@]} -ge 2 ]; then
	#### To reorder bigwig files and sample names by ChIPseq mark
	sorted_marks=()
	sorted_labels=()
	for mark in ${uniq_chip_mark_list[@]}
	do
		for sample in ${chip_sample_list[@]}
		do
			if [[ ${sample} =~ ${mark} ]]; then
				sorted_labels+=("${sample}")
			fi
		done
		for bw in ${chip_bw_list[@]}
		do
			if [[ ${bw} =~ ${mark} ]]; then
				sorted_marks+=("${bw}")
			fi
		done
	done
	if [ ${#sorted_marks[@]} -ge 1 ]; then
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
				awk -v OFS="\t" -v d=${DEG} '$11==d' ${file} > combined/matrix/temp_regions_${analysisname}_DEG_${filename}_${DEG}.bed
			done
			regions_files+=("combined/matrix/temp_regions_${analysisname}_DEG_${filename}_UP.bed" "combined/matrix/temp_regions_${analysisname}_DEG_${filename}_DOWN.bed")
			nbup=$(wc -l combined/matrix/temp_regions_${analysisname}_DEG_${filename}_UP.bed | awk '{print $1}')
			nbdown=$(wc -l combined/matrix/temp_regions_${analysisname}_DEG_${filename}_DOWN.bed | awk '{print $1}')
			regions_labels+=("${filename}_UP(${nbup})" "${filename}_DOWN(${nbdown})")
			printf "%s\t:\t%s\n" "${file}" "${regions_labels[*]}"
			awk -v OFS="\t" 'NR>1 {print $1,$2,$3}' ${file} >> combined/matrix/temp_regions_${analysisname}_all_DEGs.bed		
		done
		sort -k1,1n -k2,2n combined/matrix/temp_regions_${analysisname}_all_DEGs.bed -u > combined/matrix/temp_regions_${analysisname}_all_DEGs_unique.bed
		printf "\nComputing matrix for DEG for each sample pairs from ${analysisname}\n"
		computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${regions_files[@]} -S ${sorted_marks[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/${analysisname}_DEG.gz --quiet
		printf "\nComputing matrix for all DEGs from ${analysisname}\n"
		computeMatrix scale-regions --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${analysisname}_all_DEGs_unique.bed -S ${sorted_marks[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/${analysisname}_all_DEGs.gz --quiet
		for mark in ${uniq_chip_mark_list[@]}
		do
			selected_samples=()
			selected_labels=()
			for sample in ${sorted_labels[@]}
			do
				if [[ ${sample} =~ ${mark} ]]; then
					selected_samples+=("${sample}_merged")
					selected_labels+=("${sample}")
				fi
			done		
			computeMatrixOperations subset -m combined/matrix/${analysisname}_DEG.gz -o combined/matrix/${analysisname}_DEG_${mark}.gz --samples ${selected_samples[@]}
			printf "\nGetting scales for the DEG matrix of ${analysisname}\n"
			printf "\nPlotting ${mark} profiles for DEG for each sample pairs from ${analysisname}\n"
			printf "%s\n%s\n" "${selected_labels[*]}" "${regions_labels[*]}"
			plotProfile -m combined/matrix/${analysisname}_DEG_${mark}.gz -out combined/plots/${analysisname}_profile_DEG_${mark}.pdf --plotType 'lines' --averageType 'median' --samplesLabel ${selected_labels[@]} --regionsLabel ${regions_labels[@]} --perGroup --numPlotsPerRow 2
		done
		
		#### To plot tissue-specific DEGs
		if [ ${#rnaseq_name_list[@]} -ge 3 ]; then
			for namei in ${rnaseq_name_list[@]}
			do
				filenames="combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed"
				nbup=$(wc -l combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed | awk '{print $1}')
				nbdown=$(wc -l combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed | awk '{print $1}')
				printf "\nComputing matrix for ${namei} specific DEG from ${analysisname}\n"
				computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${filenames} -S ${sorted_marks[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o combined/matrix/${analysisname}_only_${namei}_DEG.gz --quiet
				printf "\nGetting scales for the ${tissue} specific DEG matrix of ${analysisname}\n"
				computeMatrixOperations dataRange -m combined/matrix/${analysisname}_only_${namei}_DEG.gz > combined/matrix/values_${analysisname}_only_${namei}_DEG.gz
				mins=()
				maxs=()
				for sample in ${sorted_labels[@]}
				do
					mini=$(grep ${sample} combined/matrix/values_${analysisname}_only_${namei}_DEG.gz | awk '{print $5}')
					mins+=("${mini}")
					maxi=$(grep ${sample} combined/matrix/values_${analysisname}_only_${namei}_DEG.gz | awk '{print $6}')
					maxs+=("${maxi}")
				done
				plotProfile -m combined/matrix/${analysisname}_only_${namei}_DEG.gz -out combined/plots/regions_${analysisname}_profile_only_${namei}_DEG.pdf --plotType 'lines' --averageType 'mean' --samplesLabel ${sorted_labels[@]} --regionsLabel "${namei}_UP" "${namei}_DOWN" --outFileNameData combined/matrix/values_${analysisname}_profile_only_${namei}_DEG.txt
				ymins=()
				ymaxs=()
				for sample in ${sorted_labels[@]}
				do
			 		ymini=$(grep ${sample} combined/matrix/values_${analysisname}_profile_only_${namei}_DEG.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
					ymaxi=$(grep ${sample} combined/matrix/values_${analysisname}_profile_only_${namei}_DEG.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m*1.2}')
					ymins+=("${ymini}")
					ymaxs+=("${ymaxi}")
				done
				printf "\nPlotting heatmap for ${namei} specific DEG from ${analysisname}\n"
				plotHeatmap -m combined/matrix/${analysisname}_only_${namei}_DEG.gz -out combined/plots/${analysisname}_heatmap_only_${namei}_DEG.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} --regionsLabel "${namei}_UP(${nbup})" "${namei}_DOWN(${nbdown})" --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --colorMap 'seismic' --interpolationMethod 'bilinear'
				plotProfile -m combined/matrix/${analysisname}_only_${namei}_DEG.gz -out combined/plots/regions_${analysisname}_profile_only_${namei}_DEG.pdf --plotType 'lines' --averageType 'mean' --samplesLabel ${sorted_labels[@]} --regionsLabel "${namei}_UP(${nbup})" "${namei}_DOWN(${nbdown})" --yMin ${ymins[@]} --yMax ${ymaxs[@]}
			done
		fi
	fi
fi
rm -f combined/matrix/*${analysisname}*

#########################################################################################
####################################### PART10 ##########################################
################# Making heatmaps on genes split by expression levels  ##################
#########################################################################################

#### To make heatmap and profile with deeptools for each tissue based on grouped expression levels (if both RNA and ChIP samples are present in a tissue)

for tissue in ${rnaseq_tissue_list[*]}
do
	if [[ " ${chip_tissue_list[@]} " =~ " ${tissue} " ]]; then
		tissue_chip_bw_list=()
		for bw in ${chip_bw_list[*]}
		do
			if [[ ${bw} =~ ${tissue} ]]; then
				tissue_chip_bw_list+=("${bw}")
			fi
		done
		tissue_labels_chip=()
		for sample in ${chip_sample_list[*]}
		do
			if [[ ${sample} =~ ${tissue} ]]; then
				tissue_labels_chip+=("${sample}")
			fi
		done
		tissue_rnaseq_bw_list_plus=()
		for bw in ${rnaseq_bw_list_plus[*]}
		do					
			if [[ ${bw} =~ ${tissue} ]]; then
				tissue_rnaseq_bw_list_plus+=("${bw}")
			fi
		done
		tissue_rnaseq_bw_list_minus=()
		for bw in ${rnaseq_bw_list_minus[*]}
		do					
			if [[ ${bw} =~ ${tissue} ]]; then
				tissue_rnaseq_bw_list_minus+=("${bw}")
			fi
		done
		tissue_labels_rna=()
		for sample in ${rnaseq_sample_list[*]}
		do
			if [[ ${sample} =~ ${tissue} ]]; then
				tissue_labels_rna+=("${sample}")
			fi
		done
		tissue_rampage_bw_list_plus=()
		for bw in ${rampage_bw_list_plus[*]}
		do					
			if [[ ${bw} =~ ${tissue} ]]; then
				tissue_rampage_bw_list_plus+=("${bw}")
			fi
		done
		tissue_rampage_bw_list_minus=()
		for bw in ${rampage_bw_list_minus[*]}
		do					
			if [[ ${bw} =~ ${tissue} ]]; then
				tissue_rampage_bw_list_minus+=("${bw}")
			fi
		done
		tissue_labels_rampage=()
		for sample in ${rampage_sample_list[*]}
		do
			if [[ ${sample} =~ ${tissue} ]]; then
				tissue_labels_rampage+=("${sample}")
			fi
		done
		tissue_shrna_bw_list_plus=()
		for bw in ${shrna_bw_list_plus[*]}
		do					
			if [[ ${bw} =~ ${tissue} ]]; then
				tissue_shrna_bw_list_plus+=("${bw}")
			fi
		done
		tissue_shrna_bw_list_minus=()
		for bw in ${shrna_bw_list_minus[*]}
		do					
			if [[ ${bw} =~ ${tissue} ]]; then
				tissue_shrna_bw_list_minus+=("${bw}")
			fi
		done
		tissue_labels_shrna=()
		for sample in ${shrna_sample_list[*]}
		do
			if [[ ${sample} =~ ${tissue} ]]; then
				tissue_labels_shrna+=("${sample}")
			fi
		done
		printf "Clustering genes by expression levels for ${tissue}\n"
		cols=($(awk -v ORS=" " -v t=${tissue} 'NR==1 {for(i=1;i<=NF;i++) if ($i~t) print i}' combined/DEG/counts_${analysisname}.txt))
		reps=${#cols[@]}
		awk -v d="${cols}" -v t=${reps} 'BEGIN {split(d, a, " ")} NR > 1 {b=0; for (i in a) b+=$(a[i]); c=b/t; print $1,c}' combined/DEG/counts_${analysisname}.txt > combined/DEG/temp_counts_${analysisname}_${tissue}.txt
		if [ -s combined/DEG/temp_expression_${analysisname}_${tissue}.bed ]; then
			rm -f combined/DEG/temp_expression_${analysisname}_${tissue}.bed
		fi
		while read ID exp
		do
			grep "${ID}" ${regionfile} | awk -v OFS="\t" -v c=${exp} '( $1 ~ /^[0-9]/ ) || ( $1 ~ /^chr[0-9]*$/ ) || ( $1 ~ /^Chr[0-9]*$/ ) {l=$3-$2; v=1000*c/l; print $1,$2,$3,".",v,$6,$4}' >> combined/DEG/temp_expression_${analysisname}_${tissue}.bed
		done < combined/DEG/temp_counts_${analysisname}_${tissue}.txt
		if [[ ${ref} == "B73_v4" ]]; then
			sort -k5,5gr combined/DEG/temp_expression_${analysisname}_${tissue}.bed | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$8,$5,$6}' > combined/DEG/sorted_expression_${analysisname}_${tissue}.bed
		else
			sort -k5,5gr combined/DEG/temp_expression_${analysisname}_${tissue}.bed | awk -F"[=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$8,$5,$6}' > combined/DEG/sorted_expression_${analysisname}_${tissue}.bed
		fi
				
		awk -v OFS="\t" -v a=${analysisname} -v b=${tissue} '{if ($5==0) printf $0"\n" > "combined/DEG/temp_sorted_"a"_"b"_exp0.bed"; else printf $0"\n" > "combined/DEG/temp_sorted_"a"_"b"_expA.bed"}' combined/DEG/sorted_expression_${analysisname}_${tissue}.bed
		tot=$(wc -l combined/DEG/temp_sorted_${analysisname}_${tissue}_expA.bed | awk '{print $1}')
		bin=$((tot/5))
		min=0
		max=${bin}
		for (( i = 1; i <= 5; i++ ))
		do
			awk -v n=${min} -v m=${max} 'NR>=n && NR <=m' combined/DEG/temp_sorted_${analysisname}_${tissue}_expA.bed > combined/DEG/temp_sorted_${analysisname}_${tissue}_exp${i}.bed
			min=$((min+bin))
			max=$((max+bin))
		done

		sorted_regions=()
		regions_labels=()
		for i in 1 2 3 4 5 0
		do
			case "${i}" in
				1) name="Top20%";;
				2) name="20-40%";;
				3) name="40-60%";;
				4) name="60-80%";;
				5) name="Bottom20%";;
				0) name="Not_expressed";;
			esac
			n=$(wc -l combined/DEG/temp_sorted_${analysisname}_${tissue}_exp${i}.bed | awk '{print $1}')
			regions_labels+=("${name}($n)")	
			sorted_regions+=("combined/DEG/temp_sorted_${analysisname}_${tissue}_exp${i}.bed")
			\cp -Tf combined/DEG/temp_sorted_${analysisname}_${tissue}_exp${i}.bed combined/DEG/temp_sorted_${analysisname}_${tissue}_exp${i}.txt
		done	
		for strand in plus minus
		do
			case "${strand}" in
				plus) 	bw_list="${tissue_chip_bw_list[*]} ${tissue_rnaseq_bw_list_plus[*]} ${tissue_rampage_bw_list_plus[*]} ${tissue_shrna_bw_list_plus[*]}"
					sign="+";;
				minus) 	bw_list="${tissue_chip_bw_list[*]} ${tissue_rnaseq_bw_list_minus[*]} ${tissue_rampage_bw_list_minus[*]} ${tissue_shrna_bw_list_minus[*]}"
					sign="-";;
			esac
			regions=()
			for i in 1 2 3 4 5 0
			do
				awk -v OFS="\t" -v s=$sign '$6==s' combined/DEG/temp_sorted_${analysisname}_${tissue}_exp${i}.txt > combined/DEG/temp_sorted_${analysisname}_${tissue}_exp${i}.bed
			done
			printf "\nComputing scale-regions $strand strand matrix for ${tissue} in ${analysisname}\n"
			computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/regions_${analysisname}_${strand}.gz --quiet
			printf "\nComputing reference-point on TSS $strand strand matrix for ${tissue} in $analysisname\n"
			computeMatrix reference-point --referencePoint "TSS" --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 8000 -p ${threads} -o combined/matrix/tss_${analysisname}_${strand}.gz --quiet
		done

		### Merging stranded matrix, extracting scales and plotting heatmaps
		for matrix in regions tss
		do
			tissue_labels="${tissue_labels_chip[*]} ${tissue_labels_rna[*]} ${tissue_labels_rampage[*]} ${tissue_labels_shrna[*]}"
			printf "\nMerging stranded matrices aligned by ${matrix} for ${tissue} in ${analysisname}\n"
			computeMatrixOperations rbind -m combined/matrix/${matrix}_${analysisname}_plus.gz combined/matrix/${matrix}_${analysisname}_minus.gz -o combined/matrix/${matrix}_${analysisname}.gz
			printf "\nGetting scales for ${matrix} matrix for ${tissue} in ${analysisname}\n"
			computeMatrixOperations dataRange -m combined/matrix/${matrix}_${analysisname}.gz > combined/matrix/values_${matrix}_${analysisname}.txt
			mins=()
			maxs=()
			for sample in ${tissue_labels[@]}
			do
				mini=$(grep ${sample} combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $5}')
				maxi=$(grep ${sample} combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $6}')
				test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
				if [[ ${test} == "yes" ]]; then
					mins+=("0")
					maxs+=("0.01")
				else
					mins+=("${mini}")
					maxs+=("${maxi}")
				fi
			done
			computeMatrixOperations sort -m combined/matrix/${matrix}_${analysisname}.gz -R ${sorted_regions[@]} -o combined/matrix/final_${matrix}_${analysisname}.gz
			plotProfile -m combined/matrix/final_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_profile_${matrix}.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --outFileNameData combined/matrix/values_${matrix}_${tissue}_${analysisname}.txt
			ymins=()
			ymaxs=()
			for sample in ${tissue_labels[@]}
			do
			 	ymini=$(grep ${sample} combined/matrix/values_${matrix}_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
				ymaxi=$(grep ${sample} combined/matrix/values_${matrix}_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m*1.2}')
				test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
				if [[ ${test} == "yes" ]]; then
					ymins+=("0")
					ymaxs+=("0.01")
				else
					ymins+=("${ymini}")
					ymaxs+=("${ymaxi}")
				fi
			done
			printf "\nPlotting heatmap for ${matrix} matrix for ${tissue} in ${analysisname} scaling by sample\n"
			plotHeatmap -m combined/matrix/final_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_heatmap_${matrix}.pdf --sortRegions keep --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
			if [[ ${matrix} == "tss" ]]; then
				printf "\nPlotting heatmap for ${matrix} matrix for ${tissue} in ${analysisname} scaling by sample by sorting by length\n"
				plotHeatmap -m combined/matrix/final_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_heatmap_${matrix}_v2.pdf --sortRegions descend --sortUsing region_length --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
			fi
			printf "\nPlotting profile for ${matrix} matrix for ${tissue} in ${analysisname} scaling by sample\n"
			plotProfile -m combined/matrix/final_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_profile_${matrix}.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --yMin ${ymins[@]} --yMax ${ymaxs[@]}
		done
	fi
done

rm -f combined/matrix/*${analysisname}*
rm -f combined/DEG/temp*${analysisname}*

#########################################################################################
####################################### PART11 ##########################################
############ Making heatmaps on distal H3K27ac peaks split by ChIP enrichment  ##########
#########################################################################################

#### To make heatmap and profile with deeptools for each tissue based on grouped H3K27ac levels at distal elements (>2kb)

uniq_chip_tissue_list=($(printf "%s\n" "${chip_tissue_list[@]}" | sort -u))

for tissue in ${uniq_chip_tissue_list[@]}
do
	tissue_labels=()
	tissue_bw_plus=()
	tissue_bw_minus=()
	test_k27ac="no"
	for sample in ${chip_sample_list[@]}
	do
		if [[ "${sample}" =~ "${tissue}_H3K27ac" ]]; then
			test_k27ac="yes"
		fi
		if [[ ${sample} =~ ${tissue} ]]; then
			tissue_labels+=("${sample}")
		fi
	done
	for bw in ${chip_bw_list[@]}
	do
		if [[ ${bw} =~ ${tissue} ]]; then
			tissue_bw_plus+=("${bw}")
			tissue_bw_minus+=("${bw}")
		fi
	done
	for bw in ${rnaseq_bw_list_plus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			tissue_bw_plus+=("${bw}")
		fi
	done
	for bw in ${rnaseq_bw_list_minus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			tissue_bw_minus+=("${bw}")
		fi
	done
	for sample in ${rnaseq_sample_list[*]}
	do
		if [[ ${sample} =~ ${tissue} ]]; then
			tissue_labels+=("${sample}")
		fi
	done
	for bw in ${rampage_bw_list_plus[*]}
	do					
		if [[ $bw =~ $tissue ]]; then
			tissue_bw_plus+=("$bw")
		fi
	done
	for bw in ${rampage_bw_list_minus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			tissue_bw_minus+=("${bw}")
		fi
	done
	for sample in ${rampage_sample_list[*]}
	do
		if [[ ${sample} =~ ${tissue} ]]; then
			tissue_labels+=("${sample}")
		fi
	done
	for bw in ${shrna_bw_list_plus[*]}
	do					
		if [[ $bw =~ $tissue ]]; then
			tissue_bw_plus+=("$bw")
		fi
	done
	for bw in ${shrna_bw_list_minus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			tissue_bw_minus+=("${bw}")
		fi
	done
	for sample in ${shrna_sample_list[*]}
	do
		if [[ ${sample} =~ ${tissue} ]]; then
			tissue_labels+=("${sample}")
		fi
	done
	
	if [[ ${test_k27ac} == "yes" ]] && [[ ${#tissue_bw_plus[@]} -ge 2 ]]; then
		printf "\nMaking heatmaps of distal enhancers (H3K27ac peak >2kb from TSS) in tissue ${tissue}\n"
		printf "\nGetting bed file of distal enhancers for tissue ${tissue}\n"
		bedtools sort -g ${ref_dir}/chrom.sizes -i ChIP/peaks/best_peaks_${line}_${tissue}_H3K27ac.bed > combined/peaks/temp_${analysisname}_${tissue}.bed
		bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/DEG/sorted_expression_${analysisname}_${tissue}.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '($1~/^[0-9]/ || $1~/^chr[0-9]/ ) {if ($17>= 2000 && $16=="+") print $1,$2+$10,$12,$4,$5,$16,$14,$15; else if ($17<= -2000 && $16=="-") print $1,$13,$2+$10,$4,$5,$16,$14,$15}' | sort -k5,5nr > combined/peaks/distal_${analysisname}_${tissue}.bed
		tot=$(wc -l combined/peaks/distal_${analysisname}_${tissue}.bed | awk '{print $1}')
		bin=$((tot/5))
		min=0
		max=${bin}
		for (( i = 1; i <= 5; i++ ))
		do
			awk -v n=${min} -v m=${max} 'NR>=n && NR <=m' combined/peaks/distal_${analysisname}_${tissue}.bed > combined/peaks/temp_distal_${analysisname}_${tissue}_group${i}.bed
			min=$((min+bin))
			max=$((max+bin))
		done

		sorted_regions=()
		regions_labels=()
		for i in 1 2 3 4 5
		do
			case "${i}" in
				1) name="Top20%";;
				2) name="20-40%";;
				3) name="40-60%";;
				4) name="60-80%";;
				5) name="Bottom20%";;
			esac
			n=$(wc -l combined/peaks/temp_distal_${analysisname}_${tissue}_group${i}.bed | awk '{print $1}')
			regions_labels+=("${name}($n)")	
			sorted_regions+=("combined/peaks/temp_distal_${analysisname}_${tissue}_group${i}.bed")
			\cp -Tf combined/peaks/temp_distal_${analysisname}_${tissue}_group${i}.bed combined/peaks/temp_distal_${analysisname}_${tissue}_group${i}.txt
		done	
		for strand in plus minus
		do
			case "${strand}" in
				plus) 	bw_list="${tissue_bw_plus[*]}"
					sign="+";;
				minus) 	bw_list="${tissue_bw_minus[*]}"
					sign="-";;
			esac
			regions=()
			for i in 1 2 3 4 5
			do
				awk -v OFS="\t" -v s=${sign} '$6==s' combined/peaks/temp_distal_${analysisname}_${tissue}_group${i}.txt > combined/peaks/temp_distal_${analysisname}_${tissue}_group${i}.bed
			done
			printf "\nComputing scale-regions ${strand} strand matrix for tissue ${tissue}\n"
			computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/regions_${analysisname}_distal_${strand}.gz --quiet
		done
		printf "\nMerging stranded matrices for tissue ${tissue}\n"
		computeMatrixOperations rbind -m combined/matrix/regions_${analysisname}_distal_plus.gz combined/matrix/regions_${analysisname}_distal_minus.gz -o combined/matrix/regions_${analysisname}_distal.gz
		printf "\nGetting scales for tissue ${tissue}\n"
		computeMatrixOperations dataRange -m combined/matrix/regions_${analysisname}_distal.gz > combined/matrix/values_regions_distal_${analysisname}.txt
		mins=()
		maxs=()
		for sample in ${tissue_labels[@]}
		do
			mini=$(grep ${sample} combined/matrix/values_regions_distal_${analysisname}.txt | awk '{print $5}')
			maxi=$(grep ${sample} combined/matrix/values_regions_distal_${analysisname}.txt | awk '{print $6}')
			test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
			if [[ ${test} == "yes" ]]; then
				mins+=("0")
				maxs+=("0.01")
			else
				mins+=("${mini}")
				maxs+=("${maxi}")
			fi
		done
		computeMatrixOperations sort -m combined/matrix/regions_${analysisname}_distal.gz -R ${sorted_regions[@]} -o combined/matrix/final_regions_${analysisname}_distal.gz
		plotProfile -m combined/matrix/final_regions_${analysisname}_distal.gz -out combined/plots/distal_${tissue}_${analysisname}_profile_mean.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --outFileNameData combined/matrix/values_regions_distal_${analysisname}.txt
		ymins=()
		ymaxs=()
		for sample in ${tissue_labels[@]}
		do
		 	ymini=$(grep ${sample} combined/matrix/values_regions_distal_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
			ymaxi=$(grep ${sample} combined/matrix/values_regions_distal_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m*1.2}')
			test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
			if [[ ${test} == "yes" ]]; then
				ymins+=("0")
				ymaxs+=("0.01")
			else
				ymins+=("${ymini}")
				ymaxs+=("${ymaxi}")
			fi
		done
		printf "\nPlotting heatmap for tissue ${tissue} in ${analysisname} scaling by sample\n"
		plotHeatmap -m combined/matrix/final_regions_${analysisname}_distal.gz -out combined/plots/distal_${tissue}_${analysisname}_heatmap_mean.pdf --sortRegions keep --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --startLabel "enhancer" --endLabel "TSS"
		printf "\nPlotting mean profile for tissue ${tissue} in ${analysisname} scaling by sample\n"
		plotProfile -m combined/matrix/final_regions_${analysisname}_distal.gz -out combined/plots/distal_${tissue}_${analysisname}_profile_mean.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --yMin ${ymins[@]} --yMax ${ymaxs[@]} --startLabel "enhancer" --endLabel "TSS"
		
		plotProfile -m combined/matrix/final_regions_${analysisname}_distal.gz -out combined/plots/distal_${tissue}_${analysisname}_profile_median.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType median --outFileNameData combined/matrix/values_regions_distal_${analysisname}.txt
		ymins=()
		ymaxs=()
		for sample in ${tissue_labels[@]}
		do
		 	ymini=$(grep ${sample} combined/matrix/values_regions_distal_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
			ymaxi=$(grep ${sample} combined/matrix/values_regions_distal_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m*1.2}')
			test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
			if [[ ${test} == "yes" ]]; then
				ymins+=("0")
				ymaxs+=("0.01")
			else
				ymins+=("${ymini}")
				ymaxs+=("${ymaxi}")
			fi
		done
		printf "\nPlotting median profile for tissue ${tissue} in ${analysisname} scaling by sample\n"
		plotProfile -m combined/matrix/final_regions_${analysisname}_distal.gz -out combined/plots/distal_${tissue}_${analysisname}_profile_median.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType median --yMin ${ymins[@]} --yMax ${ymaxs[@]} --startLabel "enhancer" --endLabel "TSS"
	else
		printf "\nTissue ${tissue} will not be processed (H3K27ac is present? ${test_k27ac}\tNumber of datasets ${#tissue_labels[*]}\n"
	fi
done

rm -f combined/matrix/*${analysisname}*
rm -f combined/peaks/temp*${analysisname}*

printf "\nCombined analysis script finished successfully for ${analysisname}\n"
touch combined/chkpts/analysis_${analysisname}
