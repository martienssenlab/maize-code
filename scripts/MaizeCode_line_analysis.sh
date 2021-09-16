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
##### PART12: gather RNAseq, shRNA and RAMPAGE expression at H3K27ac distal peaks and make various scatter plots
##### Under development:
##### PART5: compare TF peaks with DEGs
##### PART6: perform differential peak/TSS calling between all pairs of RAMPAGE samples
##### PART7: perform differential peak calling between all pairs of ChIP samples for each mark

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

############################################################################################
########################################## PART2 ###########################################
################# Differential gene expression analysis between tissues ####################
############################################################################################

if [ ${#rnaseq_sample_list[@]} -ge 2 ]; then
	#### This step will need to be automatized to potentially change which line/organism to have the GO terms for 
	#### It would require people to have the required files or dowload them
	if [[ ${ref} == "B73_v4" ]] || [[ ${ref} == "W22_v2" ]]; then
		if [ ! -d combined/GO ]; then
			mkdir combined/GO
		fi
		if [ ! -d combined/GO/${ref} ]; then
			mkdir combined/GO/${ref}
		fi
		if [ ! -d combined/GO/${ref}/org.Zmays.eg.db ]; then
			if [ ! -s combined/GO/${ref}/${ref}_infoGO.tab ]; then
				printf "\nCopying GO information file\n"
				cp /grid/martienssen/data_norepl/dropbox/maizecode/GO/${ref}_infoGO.tab combined/GO/${ref}/
			fi
			if [ ! -d combined/GO/${ref}/${ref}_genes_info.tab ]; then
				printf "\nCopying gene information file\n"
				cp -r /grid/martienssen/data_norepl/dropbox/maizecode/GO/${ref}_genes_info.tab combined/GO/${ref}/
			fi
			printf "\nCreating GO database\n"
			Rscript --vanilla ${mc_dir}/MaizeCode_R_build_GOdatabase.r combined/GO/${ref}/${ref}_infoGO.tab combined/GO/${ref}/${ref}_genes_info.tab ${ref}
		fi
	fi	
	#### To make a count table for all RNAseq samples in samplefile
	printf "\nPreparing count table for RNAseq samples in ${analysisname}\n"
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
				grep "gene:" RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab | sed 's/gene://' | awk -v OFS="\t" -v t=${namei} -v j=${j} 'BEGIN {print t"_Rep"j} {print $2}' > combined/DEG/col_AZ_${i}_${analysisname}_${sample}_Rep${j}.txt
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
				awk -v OFS="\t" -v t=${namei} -v j=${j} 'BEGIN {print t"_Rep"j}  $1 !~ /^N_/ {print $2}' RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab > combined/DEG/col_AZ_${i}_${analysisname}_${sample}_Rep${j}.txt
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
	if [[ ${ref} == "B73_v4" ]] || [[ ${ref} == "W22_v2" ]]; then
		printf "\nLaunching DEG analysis (and GO) with R version:\n"
		R --version
		Rscript --vanilla ${mc_dir}/MaizeCode_R_DEG_GO.r ${ref} combined/DEG/counts_${analysisname}.txt combined/DEG/samples_${analysisname}.txt ${analysisname} ${regionfile}
	else
		printf "\nLaunching DEG analysis with R version:\n"
		R --version
		Rscript --vanilla ${mc_dir}/MaizeCode_R_DEG.r combined/DEG/counts_${analysisname}.txt combined/DEG/samples_${analysisname}.txt ${analysisname} ${regionfile}
	fi
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
			printf "\nGetting DEGs specific in ${namei}\n"
			filearraydown=( $(ls combined/DEG/DEG_${analysisname}_*vs*.txt.DOWN.temp.bed | grep "${namei}") )
			filearrayup=( $(ls combined/DEG/DEG_${analysisname}_*vs*.txt.UP.temp.bed | grep "${namei}") )
			i=0
			max=$((${#filearraydown[@]}-1))
			cat ${filearraydown[i]} > combined/DEG/temp_tissue_spec_DEG_${analysisname}_DOWN_${i}.txt
			cat ${filearrayup[i]} > combined/DEG/temp_tissue_spec_DEG_${analysisname}_UP_${i}.txt
			while [ ${i} -lt ${max} ]
			do
				j=$((i+1))
				bedtools intersect -wa -a combined/DEG/temp_tissue_spec_DEG_${analysisname}_DOWN_${i}.txt -b ${filearraydown[j]} > combined/DEG/temp_tissue_spec_DEG_${analysisname}_DOWN_${j}.txt
				bedtools intersect -wa -a combined/DEG/temp_tissue_spec_DEG_${analysisname}_UP_${i}.txt -b ${filearrayup[j]} > combined/DEG/temp_tissue_spec_DEG_${analysisname}_UP_${j}.txt
				i=$((i+1))
			done
			cat combined/DEG/temp_tissue_spec_DEG_${analysisname}_DOWN_${max}.txt | sort -u > combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed
			cat combined/DEG/temp_tissue_spec_DEG_${analysisname}_UP_${max}.txt | sort -u > combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed
			rm -f combined/DEG/DEG_${analysisname}_*.temp.bed
			rm -f combined/DEG/temp_tissue_spec_DEG_${analysisname}*
			if [[ ${ref} == "B73_v4" ]] || [[ ${ref} == "W22_v2" ]]; then
				printf "\nMaking GO enrichment plot for ${namei} tissue with R version:\n"
				R --version
				Rscript --vanilla ${mc_dir}/MaizeCode_R_GO.r ${ref} combined/DEG/counts_${analysisname}.txt combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed ${namei}_DOWN_in_${analysisname}
				Rscript --vanilla ${mc_dir}/MaizeCode_R_GO.r ${ref} combined/DEG/counts_${analysisname}.txt combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed ${namei}_UP_in_${analysisname}
			fi
		done
	fi
	printf "\nCalculating DEG summary stats\n"
	numsample=${#rnaseq_name_list[@]}
	numsamplemin1=$((numsample - 1))

	if [ -e combined/reports/summary_DEG_numbers_${analysisname}.txt ]; then
		rm -f combined/reports/summary_DEG_numbers_${analysisname}.txt
	fi

	for ((i=0; i<=numsamplemin1; i++))
	do
		namei=${rnaseq_name_list[i]}
		for ((j=(i+1); j<=numsamplemin1; j++))
		do
			namej=${rnaseq_name_list[j]}
			if [ -s combined/DEG/DEG_${analysisname}_${namei}_vs_${namej}.txt ]; then
				awk -v a=${namei} -v b=${namej} -v OFS="\t" '{if ($11== "UP") c+=1; else d+=1} END {print a" vs "b,"UP:"c,"DOWN:"d}' combined/DEG/DEG_${analysisname}_${namei}_vs_${namej}.txt >> combined/reports/summary_DEG_numbers_${analysisname}.txt
			elif [ -s combined/DEG/DEG_${analysisname}_${namej}_vs_${namei}.txt ]; then
				awk -v a=${namei} -v b=${namej} -v OFS="\t" '{if ($11== "DOWN") c+=1; else d+=1} END {print a" vs "b,"UP:"c,"DOWN:"d}' combined/DEG/DEG_${analysisname}_${namej}_vs_${namei}.txt >> combined/reports/summary_DEG_numbers_${analysisname}.txt
			fi
		done
		if [ -s combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed ] && [ -s combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed ]; then
			up=$(wc -l combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed | awk '{print $1}')
			down=$(wc -l combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed | awk '{print $1}')
			awk -v a=${namei} -v u=${up} -v d=${down} -v OFS="\t" 'BEGIN {print a" only","UP:"u,"DOWN:"d}' >> combined/reports/summary_DEG_numbers_${analysisname}.txt
		elif [ -s combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed ]; then
			up=$(wc -l combined/DEG/only_${namei}_DEG_UP_${analysisname}.bed | awk '{print $1}')
			down="0"
			awk -v a=${namei} -v u=${up} -v d=${down} -v OFS="\t" 'BEGIN {print a" only","UP:"u,"DOWN:"d}' >> combined/reports/summary_DEG_numbers_${analysisname}.txt
		elif [ -s combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed ]; then
			up="0"
			down=$(wc -l combined/DEG/only_${namei}_DEG_DOWN_${analysisname}.bed | awk '{print $1}')
			awk -v a=${namei} -v u=${up} -v d=${down} -v OFS="\t" 'BEGIN {print a" only","UP:"u,"DOWN:"d}' >> combined/reports/summary_DEG_numbers_${analysisname}.txt
		fi
	done
else 
	printf "\nNo differential gene expression analysis performed (not enough RNAseq samples)\n"
fi

#############################################################################################
########################################### PART3 ###########################################
########################## Overlapping ChIPseq peaks - Upset plot  ##########################
#############################################################################################

#### To make a single file containing all overlapping peaks and plot

if [ ${#chip_sample_list[@]} -ge 1 ]; then
	printf "\nPreparing merged peaks file for $analysisname\n"
	if [ -s combined/peaks/tmp_peaks_${analysisname}.bed ]; then
		rm -f combined/peaks/tmp_peaks_${analysisname}.bed
	fi
	for sample in ${chip_sample_list[@]}
	do
		case "${sample}" in
			*H3K4me1*|*H3K27me1*) export peaktype="broad";;
			*H3K27ac*|*H3K4me3*|*H3K27me2*) export peaktype="narrow";;
		esac
		awk -v OFS="\t" -v s=${sample} '{print $1,$2,$3,s}' ChIP/peaks/selected_peaks_${sample}.${peaktype}Peak | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_${analysisname}.bed
	done
	sort -k1,1 -k2,2n combined/peaks/tmp_peaks_${analysisname}.bed > combined/peaks/tmp2_peaks_${analysisname}.bed
	bedtools merge -i combined/peaks/tmp2_peaks_${analysisname}.bed -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{print $1,$2,$3,"Peak_"NR,$4}'> combined/peaks/tmp3_peaks_${analysisname}.bed
	#### To get distance to closest gene (and the gene model name)
	printf "\nGetting closest region of $analysisname\n"
	if [[ ${ref} == "B73_v4" ]]; then
		bedtools closest -a combined/peaks/tmp3_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{print $1,$2,$3,$4,$12,".",$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/peaks/peaks_${analysisname}.bed
	else
		bedtools closest -a combined/peaks/tmp3_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{print $1,$2,$3,$4,$12,".",$5,$9}' | awk -F"[=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/peaks/peaks_${analysisname}.bed
	fi
	rm -f combined/peaks/tmp*_peaks_${analysisname}.bed
	#### To create a matrix of peak presence in each sample
	printf "\nCreating matrix file for ${analysisname}\n"
	for sample in ${chip_sample_list[@]}
	do
		printf "${sample}\n" > combined/peaks/temp_col_${analysisname}_${sample}.txt
		awk -v OFS="\t" -v s=${sample} '{if ($0 ~ s) print "1"; else print "0"}' combined/peaks/peaks_${analysisname}.bed >> combined/peaks/temp_col_${analysisname}_${sample}.txt
	done
	#### To group peaks based on their distance (gene body (x=0), promoter (0<x<2kb upstream), terminator (0<x<2kb downstream), distal)
	awk -v OFS="\t" 'BEGIN {printf "PeakID\tDistance\tGroup\n"} {if ($5<-2000) {d="Distal_downstream"; a=-$5} else if ($5<0) {d="Terminator"; a=-$5} else if ($5==0) {d="Gene_body"; a=$5} else if ($5>2000) {d="Distal_upstream"; a=$5} else {d="Promoter"; a=$5} print $4,a,d}' combined/peaks/peaks_${analysisname}.bed > combined/peaks/temp_col_${analysisname}_AAA.txt
	paste combined/peaks/temp_col_${analysisname}_*.txt | uniq > combined/peaks/matrix_upset_ChIP_${analysisname}.txt
	rm -f combined/peaks/temp_col_${analysisname}_*.txt
	#### To make an Upset plot highlighting peaks in gene bodies
	printf "\nCreating Upset plot for ${analysisname} with R version:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_Upset_ChIP.r combined/peaks/matrix_upset_ChIP_${analysisname}.txt ${analysisname}
fi

#############################################################################################
########################################### PART4 ###########################################
################## Overlapping TF peaks (w/ or w/o H3K27ac) - Upset plot  ###################
#############################################################################################

#### To make a single file containing all H3K27ac peaks of the same analysis or in the same line if possible

if [ -s combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed ]; then
	rm -f combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed
fi

nfile=0
prev_tissues=()
for file in ChIP/peaks/selected_peaks_${line}_*_H3K27ac.narrowPeak
do
	if [ -e "${file}" ]; then
		tmp1=${file##*/selected_peaks_${line}_}
		tissue=${tmp1%%_H3K27ac.narrowPeak}
		awk -v OFS="\t" -v t=${tissue} '{print $1,$2,$3,t}' ${file} | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed
		nfile=$((nfile+1))
		prev_tissues+=("${tissue}")
	fi
done

if [ ${#chip_sample_list[@]} -ge 1 ]; then
	printf "\nPreparing merged H3K27ac peaks from files in $analysisname\n"
	for sample in ${chip_sample_list[@]}
	do
		if [[ "${sample}" == *H3K27ac* ]]; then
			awk -v OFS="\t" -v s=${sample} '{print $1,$2,$3,s}' ChIP/peaks/selected_peaks_${sample}.narrowPeak | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed
		fi
	done
	sort -k1,1 -k2,2n combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed > combined/peaks/tmp2_peaks_H3K27ac_${analysisname}.bed
	bedtools merge -i combined/peaks/tmp2_peaks_H3K27ac_${analysisname}.bed | sort -k1,1 -k2,2n > combined/peaks/merged_peaks_H3K27ac_${analysisname}.bed
	rm -f combined/peaks/tmp*_peaks_H3K27ac_${analysisname}.bed
	k27file="yes"
elif [ ${nfile} -gt 0 ]; then
	printf "\nPreparing merged H3K27ac peaks from previously analyzed files, containing the following tissue(s):\n${prev_tissues[*]}"
	sort -k1,1 -k2,2n combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed > combined/peaks/tmp2_peaks_H3K27ac_${analysisname}.bed
	bedtools merge -i combined/peaks/tmp2_peaks_H3K27ac_${analysisname}.bed | sort -k1,1 -k2,2n > combined/peaks/merged_peaks_H3K27ac_${analysisname}.bed
	rm -f combined/peaks/tmp*_peaks_H3K27ac_${analysisname}.bed
	k27file="yes"
else
	printf "\nNo H3K27ac files found for ${line} line\n"
	k27file="no"
fi

#### To make a single file containing all TF peaks of the same analysis with or without H3K27ac peaks

if [ ${#tf_sample_list[@]} -ge 1 ]; then
	printf "\nPreparing merged TF peaks file for ${analysisname}\n"
	if [ -s combined/peaks/tmp_peaks_${analysisname}.bed ]; then
		rm -f combined/peaks/tmp_peaks_${analysisname}.bed
	fi
	if [[ "${k27file}" == "yes" ]]; then
		for sample in ${tf_sample_list[@]} H3K27ac
		do
			case "${sample}" in
				H3K27ac)	file="combined/peaks/merged_peaks_H3K27ac_${analysisname}.bed";;
				*)	file="TF/peaks/idr_${sample}.narrowPeak";;
			esac
			awk -v OFS="\t" -v s=${sample} '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {print $1,$2,$3,s}' ${file} | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_${analysisname}.bed
		done
	else 
		for sample in ${tf_sample_list[@]}
		do
			awk -v OFS="\t" -v s=${sample} '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {print $1,$2,$3,s}' TF/peaks/idr_${sample}.narrowPeak | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_${analysisname}.bed
		done
	fi
	sort -k1,1 -k2,2n combined/peaks/tmp_peaks_${analysisname}.bed > combined/peaks/tmp2_peaks_${analysisname}.bed
	bedtools merge -i combined/peaks/tmp2_peaks_${analysisname}.bed -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" ' $4 != "H3K27ac" {print $1,$2,$3,"Peak_"NR,$4}'> combined/peaks/tmp3_peaks_${analysisname}.bed
	#### To get distance to closest gene (and the gene model name)
	printf "\nGetting closest region of ${analysisname}\n"
	if [[ ${ref} == "B73_v4" ]]; then
		bedtools closest -a combined/peaks/tmp3_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{print $1,$2,$3,$4,$12,".",$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/peaks/peaks_${analysisname}.bed
	else
		bedtools closest -a combined/peaks/tmp3_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{print $1,$2,$3,$4,$12,".",$5,$9}' | awk -F"[:=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/peaks/peaks_${analysisname}.bed
	fi
	rm -f combined/peaks/tmp*_peaks_${analysisname}.bed
	#### To create a matrix of peak presence in each sample
	printf "\nCreating matrix file for ${analysisname}\n"
	if [[ "${k27file}" == "yes" ]]; then	
		for sample in ${tf_sample_list[@]} H3K27ac
		do
			printf "${sample}\n" > combined/peaks/temp_col_${analysisname}_${sample}.txt
			awk -v OFS="\t" -v s=${sample} '{if ($0 ~ s) print "1"; else print "0"}' combined/peaks/peaks_${analysisname}.bed >> combined/peaks/temp_col_${analysisname}_${sample}.txt
		done
	else
		for sample in ${tf_sample_list[@]}
		do
			printf "${sample}\n" > combined/peaks/temp_col_${analysisname}_${sample}.txt
			awk -v OFS="\t" -v s=${sample} '{if ($0 ~ s) print "1"; else print "0"}' combined/peaks/peaks_${analysisname}.bed >> combined/peaks/temp_col_${analysisname}_${sample}.txt
		done
	fi
	#### To group peaks based on their distance (gene body (x=0), promoter (0<x<2kb upstream), terminator (0<x<2kb downstream), distal)
	awk -v OFS="\t" 'BEGIN {printf "PeakID\tDistance\tGroup\n"} {if ($5<-2000) {d="Distal_downstream"; a=-$5} else if ($5<0) {d="Terminator"; a=-$5} else if ($5==0) {d="Gene_body"; a=$5} else if ($5>2000) {d="Distal_upstream"; a=$5} else {d="Promoter"; a=$5} print $4,a,d}' combined/peaks/peaks_${analysisname}.bed > combined/peaks/temp_col_${analysisname}_AAA.txt
	paste combined/peaks/temp_col_${analysisname}_*.txt | uniq > combined/peaks/matrix_upset_TF_${analysisname}.txt
	rm -f combined/peaks/temp_col_${analysisname}_*.txt
	#### To make an Upset plot highlighting peaks in gene bodies
	printf "\nCreating Upset plot for ${analysisname} with R version:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_Upset_TF.r ${k27file} ${analysisname} combined/peaks/matrix_upset_TF_${analysisname}.txt
fi

#############################################################################################
########################################### PART5 ###########################################
################################# TF peaks vs DEG - Bar plot  ###############################
#############################################################################################

#### To check relationship between TF peaks and DEG genes

############################################################################################
########################################## PART6 ###########################################
############################# Overlapping TSS - Upset plot  ################################
############################################################################################

#### To make a single file containing all overlapping TSS

############################################################################################
########################################## PART7 ###########################################
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
		# if [[ ${sample} =~ ${mark} ]]; then
			# sample_line_mark+=("${sample}")
		# fi
	# done
	# numsample=${#sample_line_mark[@]}
	# numsamplemin1=$((numsample-1))
	# numsamplemin2=$((numsample-2))
	# case "${mark}" in
		# H3K4me1|H3K27me2) 	peaktype="broad"
					# l=500
					# g=400;;
		# H3K27ac|H3K4me3|H3K27me2) 	peaktype="narrow"
					# l=250
					# g=150;;
	# esac
	# if [ ${numsample} -ge 2 ]; then
		# i=0
		# while [ ${i} -lt ${numsamplemin2} ]
		# do
			# a1=$(grep "properly paired" ChIP/reports/flagstat_${sample_line_mark[i]}_Rep1.txt | awk '{print $1}')
			# a2=$(grep "properly paired" ChIP/reports/flagstat_${sample_line_mark[i]}_Rep2.txt | awk '{print $1}')
			# a=$((a1+a2))
			# j=$((i+1))
			# while [ ${j} -le ${numsamplemin1} ]
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
		# printf "No differential ChIPseq peak analysis performed on ${mark} (not enough samples)\n"
	# fi
# done

#########################################################################################
####################################### PART8 ###########################################
########################## Making heatmaps on all genes  ################################
#########################################################################################

#### To make heatmaps and profiles with deeptools
#### By default, it does both scale-regions and reference-point on start of bedfile provided
#### By default, it does heatmap on all the data and corresponding profiles

if [[ "${total}" == "No" ]]; then
	printf "\nPartial combined analysis script finished successfully for ${analysisname}\n"
	touch combined/chkpts/analysis_${analysisname}
	exit 0
fi

total_sample_number=$((${#chip_sample_list[@]} + ${#rnaseq_sample_list[@]} + ${#rampage_sample_list[@]} + ${#shrna_sample_list[@]}))
if [ ${total_sample_number} -lt 2 ]; then
	printf "\nNot enough samples for deeptools analysis for ${analysisname}\nAnalysis is thus finished!\n"
	touch combined/chkpts/analysis_${analysisname}
	exit 0
fi

printf "\nDoing analysis for ${analysisname} with deeptools version:\n"
deeptools --version

uniq_chip_mark_list=($(printf "%s\n" "${chip_mark_list[@]}" | sort -u))

#### Splitting the region file by strand
awk -v OFS="\t" '$6=="+"' ${regionfile} > combined/matrix/temp_regions_${regionname}_plus.bed
awk -v OFS="\t" '$6=="-"' ${regionfile} > combined/matrix/temp_regions_${regionname}_minus.bed

#### Reordering the samples by ChIPseq mark
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

#### Computing the stranded matrix
for strand in plus minus
do
	case "${strand}" in
		plus) 	bw_list="${sorted_marks[@]} ${rnaseq_bw_list_plus[@]} ${rampage_bw_list_plus[@]} ${shrna_bw_list_plus[@]}";;
		minus) 	bw_list="${sorted_marks[@]} ${rnaseq_bw_list_minus[@]} ${rampage_bw_list_minus[@]} ${shrna_bw_list_minus[@]}";;
	esac
	printf "\nComputing scale-regions ${strand} strand matrix for ${analysisname}\n"
	computeMatrix scale-regions -q --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${regionname}_${strand}.bed -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/regions_${analysisname}_${strand}.gz
	printf "\nComputing reference-point on TSS ${strand} strand matrix for ${analysisname}\n"
	computeMatrix reference-point --referencePoint "TSS" -q --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${regionname}_${strand}.bed -S ${bw_list} -bs 50 -b 2000 -a 8000 -p ${threads} -o combined/matrix/tss_${analysisname}_${strand}.gz
done
rm -f combined/matrix/temp_regions_${regionname}_*.bed

#### Merging stranded matrix, extracting scales and plotting heatmaps
all_samples=()
all_lables=()
if [ ${#sorted_marks[@]} -gt 0 ]; then
	printf "\nIncluding ChIPseq samples\n"
	all_samples+=("${uniq_chip_mark_list[*]}")
	all_labels+=("${sorted_labels[*]}")
fi
if [ ${#rnaseq_bw_list_plus[@]} -gt 0 ]; then
	printf "\nIncluding RNAseq samples\n"
	all_samples+=("RNAseq")
	all_labels+=("${rnaseq_sample_list[*]}")
fi
if [ ${#rampage_bw_list_plus[@]} -gt 0 ]; then
	printf "\nIncluding RAMPAGE samples\n"
	all_samples+=("RAMPAGE")
	all_labels+=("${rampage_sample_list[*]}")
fi
if [ ${#shrna_bw_list_plus[@]} -gt 0 ]; then
	printf "\nIncluding shRNA samples\n"
	all_samples+=("shRNA")
	all_labels+=("${shrna_sample_list[*]}")
fi
for matrix in regions tss
do
	printf "\nMerging stranded matrices aligned by ${matrix} of ${analysisname}\n"
	computeMatrixOperations rbind -m combined/matrix/${matrix}_${analysisname}_plus.gz combined/matrix/${matrix}_${analysisname}_minus.gz -o combined/matrix/${matrix}_${analysisname}.gz
	printf "\nGetting scales for ${matrix} matrix of ${analysisname}\n"
	computeMatrixOperations dataRange -m combined/matrix/${matrix}_${analysisname}.gz > combined/matrix/values_${matrix}_${analysisname}.txt
	plotProfile -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/temp_${matrix}_${analysisname}_profile.pdf --samplesLabel ${all_labels[@]} --averageType mean --outFileNameData combined/matrix/values_profile_${matrix}_${analysisname}.txt
	rm -f combined/plots/temp_${matrix}_${analysisname}_profile.pdf
	mins=()
	maxs=()
	ymins=()
	ymaxs=()
	for mark in ${all_samples[@]}
	do
		mini=$(grep "${mark}" combined/matrix/values_${matrix}_${analysisname}.txt | awk 'BEGIN {m=999999} {a=$5; if (a<m) m=a;} END {print m}')
		maxi=$(grep "${mark}" combined/matrix/values_${matrix}_${analysisname}.txt | awk 'BEGIN {m=-999999} {a=$6; if (a>m) m=a;} END {print m}')
		num=$(grep "${mark}" combined/matrix/values_${matrix}_${analysisname}.txt | wc -l)
		test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		if [[ ${test} == "yes" ]]; then
			mini=("0")
			maxi=("0.01")
		fi
		for i in $(seq 1 ${num})
		do
			mins+=("${mini}")
			maxs+=("${maxi}")
		done		
		ymini=$(grep "${mark}" combined/matrix/values_profile_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
		ymaxi=$(grep "${mark}" combined/matrix/values_profile_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m*1.2}')
		num=$(grep "${mark}" combined/matrix/values_profile_${matrix}_${analysisname}.txt | wc -l)
		test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		if [[ ${test} == "yes" ]]; then
			ymini=("0")
			ymaxi=("0.01")
		fi
		for i in $(seq 1 ${num})
		do
			ymins+=("${ymini}")
			ymaxs+=("${ymaxi}")
		done
	done
	
	mins2=()
	maxs2=()
	for sample in ${all_labels[@]}
	do
		mini=$(grep ${sample} combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $5}')
		maxi=$(grep ${sample} combined/matrix/values_${matrix}_${analysisname}.txt | awk '{print $6}')
		test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		if [[ ${test} == "yes" ]]; then
			mins2+=("0")
			maxs2+=("0.01")
		else
			mins2+=("${mini}")
			maxs2+=("${maxi}")
		fi
	done
	ymins2=()
	ymaxs2=()
	for sample in ${all_labels[@]}
	do
		ymini=$(grep ${sample} combined/matrix/values_profile_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
		ymaxi=$(grep ${sample} combined/matrix/values_profile_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m*1.2}')
		test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		if [[ ${test} == "yes" ]]; then
			ymins2+=("0")
			ymaxs2+=("0.01")
		else
			ymins2+=("${ymini}")
			ymaxs2+=("${ymaxi}")
		fi
	done
	printf "\nPlotting heatmap for ${matrix} matrix of ${analysisname} scaling by mark\n"
	plotHeatmap -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_${matrix}.pdf --sortRegions descend --sortUsing mean --samplesLabel ${all_labels[@]} --regionsLabel ${regionname} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
	printf "\nPlotting heatmap for ${matrix} matrix of ${analysisname} scaling by sample\n"
	plotHeatmap -m combined/matrix/${matrix}_${analysisname}.gz -out combined/plots/${analysisname}_heatmap_${matrix}_v2.pdf --sortRegions descend --sortUsing mean --samplesLabel ${all_labels[@]} --regionsLabel ${regionname} --colorMap 'seismic' --zMin ${mins2[@]} --zMax ${maxs2[@]} --yMin ${ymins2[@]} --yMax ${ymaxs2[@]} --interpolationMethod 'bilinear'
done

rm -f combined/matrix/*${analysisname}*.gz

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
			printf "\nComputing scale-regions $strand strand matrix for split expression in ${tissue} in ${analysisname}\n"
			computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/regions_${analysisname}_${strand}.gz --quiet
			printf "\nComputing reference-point on TSS $strand strand matrix for split expression in ${tissue} in $analysisname\n"
			computeMatrix reference-point --referencePoint "TSS" --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 8000 -p ${threads} -o combined/matrix/tss_${analysisname}_${strand}.gz --quiet
		done

		### Merging stranded matrix, extracting scales and plotting heatmaps
		for matrix in regions tss
		do
			tissue_labels="${tissue_labels_chip[*]} ${tissue_labels_rna[*]} ${tissue_labels_rampage[*]} ${tissue_labels_shrna[*]}"
			printf "\nMerging stranded matrices aligned by ${matrix} for split expression in ${tissue} in ${analysisname}\n"
			computeMatrixOperations rbind -m combined/matrix/${matrix}_${analysisname}_plus.gz combined/matrix/${matrix}_${analysisname}_minus.gz -o combined/matrix/${matrix}_${analysisname}.gz
			printf "\nGetting scales for ${matrix} matrix for split expression in ${tissue} in ${analysisname}\n"
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
			printf "\nPlotting heatmap for ${matrix} matrix for split expression in ${tissue} in ${analysisname} scaling by sample\n"
			plotHeatmap -m combined/matrix/final_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_heatmap_${matrix}.pdf --sortRegions keep --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
			if [[ ${matrix} == "tss" ]]; then
				printf "\nPlotting heatmap for ${matrix} matrix for split expression in ${tissue} in ${analysisname} scaling by sample by sorting by length\n"
				plotHeatmap -m combined/matrix/final_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_heatmap_${matrix}_v2.pdf --sortRegions descend --sortUsing region_length --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
			fi
			printf "\nPlotting profile for ${matrix} matrix for split expression in ${tissue} in ${analysisname} scaling by sample\n"
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

h3k27actissues=()
for tissue in ${uniq_chip_tissue_list[@]}
do
	tissue_labels=()
	tissue_bw_plus=()
	tissue_bw_minus=()
	test_k27ac="no"
	for sample in ${chip_sample_list[@]}
	do
		if [[ "${sample}" =~ "${tissue}_H3K27ac" ]]; then
			h3k27actissues+=("${tissue}")
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
		bedtools sort -g ${ref_dir}/chrom.sizes -i combined/DEG/sorted_expression_${analysisname}_${tissue}.bed > combined/peaks/temp_${analysisname}_${tissue}_expression.bed
		bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '($1~/^[0-9]/ || $1~/^chr[0-9]/ ) {if ($17>= 2000 && $16=="+") print $1,$2+$10,$12,$4,$5,$16,$14,$15; else if ($17<= -2000 && $16=="-") print $1,$13,$2+$10,$4,$5,$16,$14,$15}' | sort -k5,5nr > combined/peaks/distal_${analysisname}_${tissue}.bed
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
			awk -v OFS="\t" -v g=${name} -v t=${tissue} '{print t,$4,$5,$7,$8,g }' combined/peaks/temp_distal_${analysisname}_${tissue}_group${i}.bed >> combined/peaks/temp2_distal_${analysisname}_${tissue}.txt
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
		printf "\nMerging stranded matrices for distal peaks in tissue ${tissue}\n"
		computeMatrixOperations rbind -m combined/matrix/regions_${analysisname}_distal_plus.gz combined/matrix/regions_${analysisname}_distal_minus.gz -o combined/matrix/regions_${analysisname}_distal.gz
		printf "\nGetting scales for distal peaks in tissue ${tissue}\n"
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
		printf "\nPlotting heatmap for distal peaks in tissue ${tissue} in ${analysisname} scaling by sample\n"
		plotHeatmap -m combined/matrix/final_regions_${analysisname}_distal.gz -out combined/plots/distal_${tissue}_${analysisname}_heatmap_mean.pdf --sortRegions keep --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --startLabel "enhancer" --endLabel "TSS"
		printf "\nPlotting mean profile for distal peaks in tissue ${tissue} in ${analysisname} scaling by sample\n"
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
		printf "\nPlotting median profile for distal peaks in tissue ${tissue} in ${analysisname} scaling by sample\n"
		plotProfile -m combined/matrix/final_regions_${analysisname}_distal.gz -out combined/plots/distal_${tissue}_${analysisname}_profile_median.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType median --yMin ${ymins[@]} --yMax ${ymaxs[@]} --startLabel "enhancer" --endLabel "TSS"		
	else
		printf "\nTissue ${tissue} will not be processed (H3K27ac is present? ${test_k27ac}\tNumber of datasets in ${tissue}: ${#tissue_labels[*]})\n"
	fi
done
if [[ ${test_k27ac} == "yes" ]] && [[ ${#tissue_bw_plus[@]} -ge 2 ]]; then
	printf "Tissue\tPeak_ID\tPeakQuality\tGID\tRPKM\tGroup\n" > combined/peaks/all_grouped_distal_peaks_${analysisname}_*.txt
	cat combined/peaks/temp2_distal_${analysisname}_*.txt > combined/peaks/all_grouped_distal_peaks_${analysisname}_*.txt
fi

rm -f combined/matrix/*${analysisname}*
rm -f combined/peaks/temp*${analysisname}*

#########################################################################################
####################################### PART12 ##########################################
############ Making scatter plots of RNA expression on distal H3K27ac peaks #############
#########################################################################################

#### To make scatter plots with R for each tissue based on RNAseq/RAMPAGE/shRNA and quality of H3K27ac distal peaks

for tissue in ${h3k27actissues[@]}
do
	header="Chr\tStart\tStop\tPeakID"
	awk -v OFS="\t" '{print $1,$2,$3,$4}' ChIP/peaks/best_peaks_${line}_${tissue}_H3K27ac.bed > combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt
	rnaseq=0
	rampage=0
	shrna=0
	for bw in ${rnaseq_bw_list_plus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			printf "\nGetting RNAseq plus strand coverage on H3K27ac peaks for ${tissue}\n"
			bigWigToBedGraph ${bw} combined/peaks/temp_RNAseq_plus_${analysisname}_${line}.bedGraph
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_RNAseq_plus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=$5/($3-$2); print l}' > combined/peaks/col_B_${analysisname}_${line}_${tissue}.txt
			rm -f combined/peaks/temp_RNAseq_plus_${analysisname}_${line}.bedGraph
			header="${header}\tRNAseq_plus"
			rnaseq=1
		fi
	done
	for bw in ${rnaseq_bw_list_minus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			printf "\nGetting RNAseq minus strand coverage on H3K27ac peaks for ${tissue}\n"
			bigWigToBedGraph ${bw} combined/peaks/temp_RNAseq_minus_${analysisname}_${line}.bedGraph
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_RNAseq_minus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=$5/($3-$2); print l}' > combined/peaks/col_C_${analysisname}_${line}_${tissue}.txt
			rm -f combined/peaks/temp_RNAseq_minus_${analysisname}_${line}.bedGraph
			header="${header}\tRNAseq_minus"
		fi
	done
	for bw in ${rampage_bw_list_plus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			printf "\nGetting RAMPAGE plus strand coverage on H3K27ac peaks for ${tissue}\n"
			bigWigToBedGraph ${bw} combined/peaks/temp_RAMPAGE_plus_${analysisname}_${line}.bedGraph
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_RAMPAGE_plus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=$5/($3-$2); print l}' > combined/peaks/col_D_${analysisname}_${line}_${tissue}.txt
			rm -f combined/peaks/temp_RAMPAGE_plus_${analysisname}_${line}.bedGraph
			header="${header}\tRAMPAGE_plus"
			rampage=1
		fi
	done
	for bw in ${rampage_bw_list_minus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			printf "\nGetting RAMPAGE minus strand coverage on H3K27ac peaks for ${tissue}\n"
			bigWigToBedGraph ${bw} combined/peaks/temp_RAMPAGE_minus_${analysisname}_${line}.bedGraph
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_RAMPAGE_minus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=$5/($3-$2); print l}' > combined/peaks/col_E_${analysisname}_${line}_${tissue}.txt
			rm -f combined/peaks/temp_RAMPAGE_minus_${analysisname}_${line}.bedGraph
			header="${header}\tRAMPAGE_minus"
		fi
	done
	for bw in ${shrna_bw_list_plus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			printf "\nGetting shRNA plus strand coverage on H3K27ac peaks for ${tissue}\n"
			bigWigToBedGraph ${bw} combined/peaks/temp_shRNA_plus_${analysisname}_${line}.bedGraph
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_shRNA_plus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=$5/($3-$2); print l}' > combined/peaks/col_F_${analysisname}_${line}_${tissue}.txt
			rm -f combined/peaks/temp_shrna_plus_${analysisname}_${line}.bedGraph
			header="${header}\tshRNA_plus"
			shrna=1
		fi
	done
	for bw in ${shrna_bw_list_minus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			printf "\nGetting shRNA minus strand coverage on H3K27ac peaks for ${tissue}\n"
			bigWigToBedGraph ${bw} combined/peaks/temp_shRNA_minus_${analysisname}_${line}.bedGraph
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_shRNA_minus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=$5/($3-$2); print l}' > combined/peaks/col_G_${analysisname}_${line}_${tissue}.txt
			rm -f combined/peaks/temp_shRNA_minus_${analysisname}_${line}.bedGraph
			header="${header}\tshRNA_minus"
		fi
	done	
	printf "${header}\n" > combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt
	paste combined/peaks/col_*_${analysisname}_${line}_${tissue}.txt >> combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt
	
	if [[ ${rnaseq} == 1 ]] && [[ ${rampage} == 1 ]] && [[ ${shrna} == 1 ]]; then
		included_samples="RNAseq_and_RAMPAGE_and_shRNA"
	elif [[ ${rnaseq} == 1 ]] && [[ ${rampage} == 1 ]]; then
		included_samples="RNAseq_and_RAMPAGE"
	elif [[ ${rnaseq} == 1 ]] && [[ ${shrna} == 1 ]]; then
		included_samples="RNAseq_and_shRNA"
	elif [[ ${rampage} == 1 ]] && [[ ${shrna} == 1 ]]; then
		included_samples="RAMPAGE_and_shRNA"
	elif [[ ${rnaseq} == 1 ]]; then
		included_samples="RNAseq"
	elif [[ ${rampage} == 1 ]]; then
		included_samples="RAMPAGE"
	elif [[ ${shRNA} == 1 ]]; then
		included_samples="shRNA"
	fi

	#### To plot correlation between RNA expression datasets at distal peaks
	printf "\nCreating scatter plot for ${analysisname} ${tissue} with R version:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_scatter_distal_peaks.r ${analysisname} ${tissue} ${line} ${included_samples} combined/peaks/all_grouped_distal_peaks_${analysisname}.txt combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt
done

printf "\nCombined analysis script finished successfully for ${analysisname}\n"
touch combined/chkpts/analysis_${analysisname}
