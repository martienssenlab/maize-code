#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 32
#$ -l m_mem_free=4G
#$ -l tmp_free=8G
#$ -o lineanalysis.log
#$ -j y
#$ -N lineanalysis

usage="
##### Script for Maize code inbred line data analysis, used by script MaizeCode_analyis.sh if -s was not set and region file exists
#####
##### sh MaiCode_line_analysis.sh -f samplefile -r regionfile [-t] [-z]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Sample, PE or SE, Reference genome directory
##### 	-r: bedfile containing the regions that are to be ploted over
#####	-t: If set, partial analysis will be performed (no heatmap with deeptools)
#####	-z: If set, partial analysis will be performed for testing
#####	-x: If set, heatmaps on repeats will be performed (can take a very long time)
##### 	-h: help, returns usage
##### 
##### After preparing all the data in the samplefile (PART1) it will, depending on the types of data present:
##### PART2: perform differential gene expression analysis among all RNAseq samples, potentially with gene ontology analysis (for B73_v4 and W22 for now)
##### PART3: produce an Upset plot of the intersection between all ChIP samples, highlighting peaks in the input regions
##### PART4: produce an Upset plot of the intersection between all TF samples, highlighting peaks in the H3K27ac peaks (if possible)
##### PART8: plot different heatmaps and metaplots of all ChIP and RNA samples in the samplefile on all input regions (and on TEs for B73_v4, need a TE gff file)
##### PART9: plot heatmaps and metaplots on DEGs
##### PART10: plot heatmaps and metaplots on genes split by expression levels
##### PART11: plot heatmaps and metaplots on H3K27ac distal peaks split by H3K27ac ChIPseq quality
##### PART12: gather RNAseq, shRNA and RAMPAGE expression at H3K27ac distal peaks and make various scatter plots
##### PART13: plot heatmaps and metaplots on enhancers (H3K27ac peaks) split by distance to nearest gene, sorted by RNA expression at enhancers (uses part2, part3, part11 and part12)
##### PART14: plot heatmaps and metaplots on enhancers (H3K27ac peaks) split by distance to nearest gene, sorted by RNA expression at enhancers (uses part13)
##### PART15: plot the distribution of RAMPAGE "TSS" peaks in genes and TEs for each tissue and an upset plot among all tissues (need a TE gff file)
##### PART16: plot the distribution of shRNA clusters in genes and TEs for each tissue and an upset plot among all tissues (need a TE gff file)
##### Under development:
##### PART5: compare TF peaks with DEGs
##### PART7: perform differential peak calling between all pairs of ChIP samples for each mark
##### PART17: plot heatmaps on all TEs

#####
##### Requirements: bedtools, deeptools, macs2, R (+R packages: ggplot2,limma,edgeR,dplyr,tidyr,stringr,gplots,cowplot,ComplexUpset,purrr)
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

while getopts ":f:r:tzh" opt; do
	case ${opt} in
		h) 	printf "${usage}\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		r)	export regionfile=${OPTARG};;
		t)	printf "\nPartial analysis to be performed (no heatmaps)\n" 
			export total="NO";;
		z)	printf "\nPartial analysis to be performed for testing\n" 
			export total="TEST";;
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
		ChIP*) 	datatype="ChIP"
			name=${line}_${tissue}_${sample};;
		RNAseq) datatype="RNA"
			name=${line}_${tissue}_${sample};;
		RAMPAGE) datatype="RNA"
			name=${line}_${tissue}_${sample};;
		shRNA) datatype="shRNA"
			name=${line}_${tissue}_${sample};;
		TF_*) datatype="TF"
			name=${data##TF_};;
	esac

	ref=${ref_dir##*/}
	if [[ ! "${ref_list[@]}" =~ "${ref}" ]]; then
		ref_list+=("${ref}")
		ref_dir_list+=("${ref_dir}")
	fi
	if [[ ! "${line_list[@]}" =~ "${line}" ]]; then
		line_list+=("${line}")
	fi
	if [[ "${datatype}" == "ChIP"* ]]; then
		if [ -e ${datatype}/tracks/${name}_merged.bw ]; then
			chip_bw_list+=("${datatype}/tracks/${name}_merged.bw")
		else
			chip_bw_list+=("${datatype}/tracks/${name}_Rep1.bw")
		fi
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
	if [[ ${ref} == "B73_v4" ]] || [[ ${ref} == "B73_v5" ]] || [[ ${ref} == "W22_v2" ]] || [[ ${ref} == "NC350_NAM" ]]; then
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
			if [ ! -s combined/GO/${ref}/${ref}_genes_info.tab ]; then
				printf "\nCopying gene information file\n"
				cp /grid/martienssen/data_norepl/dropbox/maizecode/GO/${ref}_genes_info.tab combined/GO/${ref}/
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
	if [[ ${ref} == "B73_v4" ]] || [[ ${ref} == "B73_v5" ]] || [[ ${ref} == "W22_v2" ]] || [[ ${ref} == "NC350_NAM" ]]; then
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
			if [ -e combined/peaks/temp_all_${namei}_${analysisname}_DEG_GID.txt ]; then
				rm -f combined/peaks/temp_all_${namei}_${analysisname}_DEG_GID.txt
			fi
			if ls combined/DEG/DEG_${analysisname}_${namei}_vs_*.txt 1> /dev/null 2>&1; then
				cat combined/DEG/DEG_${analysisname}_${namei}_vs_*.txt | awk '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/) {print $4}' > combined/peaks/temp_all_${namei}_${analysisname}_DEG_GID.txt
				for file in $(ls combined/DEG/DEG_${analysisname}_${namei}_vs_*.txt)
				do
					awk -v OFS="\t" '$11 == "DOWN" {print $1,$2,$3,$4,".",$6}' ${file} > ${file}.DOWN.temp.bed
					awk -v OFS="\t" '$11 == "UP" {print $1,$2,$3,$4,".",$6}' ${file} > ${file}.UP.temp.bed
				done
			fi
			if ls combined/DEG/DEG_${analysisname}_*_vs_${namei}.txt 1> /dev/null 2>&1; then
				cat combined/DEG/DEG_${analysisname}_*_vs_${namei}.txt | awk '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/) {print $4}' >> combined/peaks/temp_all_${namei}_${analysisname}_DEG_GID.txt
				for file in $(ls combined/DEG/DEG_${analysisname}_*_vs_${namei}.txt)
				do
					awk -v OFS="\t" '$11 == "DOWN" {print $1,$2,$3,$4,".",$6}' ${file} >> ${file}.UP.temp.bed
					awk -v OFS="\t" '$11 == "UP" {print $1,$2,$3,$4,".",$6}' ${file} >> ${file}.DOWN.temp.bed
				done
			fi
			printf "\nGetting DEGs specific in ${namei}\n"
			sort -u combined/peaks/temp_all_${namei}_${analysisname}_DEG_GID.txt > combined/peaks/all_${namei}_${analysisname}_DEG_GID.txt
			filearraydown=( $(ls combined/DEG/DEG_${analysisname}_*_vs_*.txt.DOWN.temp.bed | grep "${namei}") )
			filearrayup=( $(ls combined/DEG/DEG_${analysisname}_*_vs_*.txt.UP.temp.bed | grep "${namei}") )
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
			if [[ ${ref} == "B73_v4" ]] || [[ ${ref} == "B73_v5" ]] || [[ ${ref} == "W22_v2" ]] || [[ ${ref} == "NC350_NAM" ]]; then
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
elif [ ${#rnaseq_sample_list[@]} -eq 1 ]; then
	namei=${rnaseq_name_list[0]}
	sample=${rnaseq_sample_list[0]}
	printf "\nPreparing count table for the RNAseq samples of the one tissue (${sample}) in ${analysisname}\n"
	printf "Replicate\tSample\tColor\n" > combined/DEG/samples_${analysisname}.txt
	numreps=$(ls -1f RNA/mapped/map_${sample}_Rep*_ReadsPerGene.out.tab | wc -l)
	for ((j=1;j<=numreps;j++))
	do
		printf "${namei}_Rep${j}\t${namei}\t1\n" >> combined/DEG/samples_${analysisname}.txt
		if [ $(grep "gene:" RNA/mapped/map_${sample}_Rep1_ReadsPerGene.out.tab | wc -l) -gt 0 ]; then
			grep "gene:" RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab | sed 's/gene://' | awk -v OFS="\t" -v t=${namei} -v j=${j} 'BEGIN {print t"_Rep"j} {print $2}' > combined/DEG/col_AZ_${i}_${analysisname}_${sample}_Rep${j}.txt
			if [ $j -eq 1 ]; then
				grep "gene:" RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab | sed 's/gene://' | awk -v OFS="\t" 'BEGIN {print "gene_ID"} {print $1}' > combined/DEG/col_AA_0_${analysisname}.txt
			fi
		else
			awk -v OFS="\t" -v t=${namei} -v j=${j} 'BEGIN {print t"_Rep"j}  $1 !~ /^N_/ {print $2}' RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab > combined/DEG/col_AZ_${i}_${analysisname}_${sample}_Rep${j}.txt
			if [ $j -eq 1 ]; then
				awk -v OFS="\t" 'BEGIN {print "gene_ID"}  $1 !~ /^N_/ {print $1}' RNA/mapped/map_${sample}_Rep${j}_ReadsPerGene.out.tab > combined/DEG/col_AA_0_${analysisname}.txt
			fi
		fi
	done
	paste combined/DEG/col_A*_${analysisname}* > combined/DEG/counts_${analysisname}.txt
	rm -f combined/DEG/col_A*_${analysisname}*
else 
	printf "\nNo differential gene expression analysis performed (not enough RNAseq samples)\n"
fi

uniq_rnaseq_tissue_list=($(printf "%s\n" "${rnaseq_tissue_list[@]}" | sort -u))

if [[ ${#uniq_rnaseq_tissue_list[@]} -ge 1 ]]; then
	for tissue in ${uniq_rnaseq_tissue_list[@]}
	do
		printf "Gathering gene expression levels for ${tissue}\n"
		cols=($(awk -v ORS=" " -v t=${tissue} 'NR==1 {for(i=1;i<=NF;i++) if ($i~t) print i}' combined/DEG/counts_${analysisname}.txt))
		reps=${#cols[@]}
		awk -v OFS="\t" -v d="${cols}" -v t=${reps} 'BEGIN {split(d, a, " ")} NR > 1 {b=0; for (i in a) b+=$(a[i]); c=b/t; print $1,c}' combined/DEG/counts_${analysisname}.txt > combined/DEG/temp_counts_${analysisname}_${tissue}.txt
		if [ -s combined/DEG/temp_expression_${analysisname}_${tissue}.bed ]; then
			rm -f combined/DEG/temp_expression_${analysisname}_${tissue}.bed
		fi
		while read ID exp
		do
			grep "${ID}" ${regionfile} | awk -v OFS="\t" -v c=${exp} '{l=$3-$2; v=1000*c/l; print $1,$2,$3,".",v,$6,$4}' >> combined/DEG/temp_expression_${analysisname}_${tissue}.bed
		done < combined/DEG/temp_counts_${analysisname}_${tissue}.txt
		if [[ ${ref} == "B73_v4" ]]; then
			awk -F"[:;]" -v OFS="\t" '{print $1,$2}' combined/DEG/temp_expression_${analysisname}_${tissue}.bed | awk -v OFS="\t" -v t=${tissue} '{print $8,t,$5}' > combined/DEG/genes_rpkm_${analysisname}_${tissue}.txt
		else
			awk -F"[=;]" -v OFS="\t" '{print $1,$2}' combined/DEG/temp_expression_${analysisname}_${tissue}.bed | awk -v OFS="\t" -v t=${tissue} '{print $8,t,$5}' > combined/DEG/genes_rpkm_${analysisname}_${tissue}.txt
		fi
		rm -f combined/DEG/temp_*_${analysisname}_${tissue}.txt
	done
	printf "GID\tTissue\tRPKM\n" > combined/DEG/genes_rpkm_${analysisname}.txt
	cat combined/DEG/genes_rpkm_${analysisname}_*.txt >> combined/DEG/genes_rpkm_${analysisname}.txt
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
			*H3K4me1*|*H3K27me1*|*H3K27me2*|*H3K27me3*|*H3K9me1*|*H3K9me2*|*H3K9me3*) export peaktype="broad";;
			*H3K27ac*|*H3K4me3*) export peaktype="narrow";;
		esac
		awk -v OFS="\t" -v s=${sample} '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {print $1,$2,$3,s}' ChIP/peaks/selected_peaks_${sample}.${peaktype}Peak | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_${analysisname}.bed
	done
	sort -k1,1 -k2,2n combined/peaks/tmp_peaks_${analysisname}.bed > combined/peaks/tmp2_peaks_${analysisname}.bed
	bedtools merge -i combined/peaks/tmp2_peaks_${analysisname}.bed -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{print $1,$2,$3,"Peak_"NR,$4}'> combined/peaks/tmp3_peaks_${analysisname}.bed
	#### To get distance to closest gene (and the gene model name)
	printf "\nGetting closest region of $analysisname\n"
	if [[ ${ref} == "B73_v4" ]]; then
		bedtools closest -a combined/peaks/tmp3_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/peaks/peaks_${analysisname}.bed
	else
		bedtools closest -a combined/peaks/tmp3_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/peaks/peaks_${analysisname}.bed
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
	#### To make an Upset plot highlighting peaks in gene bodies (if more than one sample are present)
	if [ ${#chip_sample_list[@]} -ge 2 ]; then
		printf "\nCreating Upset plot for ${analysisname} with R version:\n"
		R --version
		Rscript --vanilla ${mc_dir}/MaizeCode_R_Upset_ChIP.r combined/peaks/matrix_upset_ChIP_${analysisname}.txt ${analysisname}
	fi
fi

#############################################################################################
########################################### PART4 ###########################################
################## Overlapping TF peaks (w/ or w/o H3K27ac) - Upset plot  ###################
#############################################################################################

#### To make a single file containing all H3K27ac peaks of the same analysis or in the same line if possible

if [ -s combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed ]; then
	rm -f combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed
fi

k27file="no"
insamplefile="no"
if [ ${#chip_sample_list[@]} -ge 1 ]; then	
	for sample in ${chip_sample_list[@]}
	do
		if [[ "${sample}" == *H3K27ac* ]]; then
			insamplefile="yes"
			awk -v OFS="\t" -v s=${sample} '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/) {print $1,$2,$3,s}' ChIP/peaks/selected_peaks_${sample}.narrowPeak | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed
		fi
	done
	if [[ "${insamplefile}" == "yes" ]]; then
		printf "\nPreparing merged H3K27ac peaks from files in ${analysisname}\n"
		sort -k1,1 -k2,2n combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed > combined/peaks/tmp2_peaks_H3K27ac_${analysisname}.bed
		bedtools merge -i combined/peaks/tmp2_peaks_H3K27ac_${analysisname}.bed | sort -k1,1 -k2,2n > combined/peaks/merged_peaks_H3K27ac_${analysisname}.bed
		rm -f combined/peaks/tmp*_peaks_H3K27ac_${analysisname}.bed
		k27file="yes"
	fi
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
if [[ "${insamplefile}" == "no" ]] && [ ${nfile} -gt 0 ]; then
	printf "\nPreparing merged H3K27ac peaks from previously analyzed files, containing the following tissue(s):\n${prev_tissues[*]}"
	sort -k1,1 -k2,2n combined/peaks/tmp_peaks_H3K27ac_${analysisname}.bed > combined/peaks/tmp2_peaks_H3K27ac_${analysisname}.bed
	bedtools merge -i combined/peaks/tmp2_peaks_H3K27ac_${analysisname}.bed | sort -k1,1 -k2,2n > combined/peaks/merged_peaks_H3K27ac_${analysisname}.bed
	rm -f combined/peaks/tmp*_peaks_H3K27ac_${analysisname}.bed
	k27file="yes"
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
				*)	file="TF/peaks/idr_${line}_${sample}.narrowPeak";;
			esac
			awk -v OFS="\t" -v s=${sample} '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {print $1,$2,$3,s}' ${file} | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_${analysisname}.bed
		done
	else 
		for sample in ${tf_sample_list[@]}
		do
			awk -v OFS="\t" -v s=${sample} '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {print $1,$2,$3,s}' TF/peaks/idr_${line}_${sample}.narrowPeak | sort -k1,1 -k2,2n -u >> combined/peaks/tmp_peaks_${analysisname}.bed
		done
	fi
	sort -k1,1 -k2,2n combined/peaks/tmp_peaks_${analysisname}.bed > combined/peaks/tmp2_peaks_${analysisname}.bed
	bedtools merge -i combined/peaks/tmp2_peaks_${analysisname}.bed -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" ' $4 != "H3K27ac" {print $1,$2,$3,"Peak_"NR,$4}'> combined/peaks/tmp3_peaks_${analysisname}.bed
	#### To get distance to closest gene (and the gene model name)
	printf "\nGetting closest region of ${analysisname}\n"
	if [[ ${ref} == "B73_v4" ]]; then
		bedtools closest -a combined/peaks/tmp3_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/peaks/TF_peaks_${analysisname}.bed
	else
		bedtools closest -a combined/peaks/tmp3_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/peaks/TF_peaks_${analysisname}.bed
	fi
	rm -f combined/peaks/tmp*_peaks_${analysisname}.bed
	#### To create a matrix of peak presence in each sample
	printf "\nCreating matrix file for ${analysisname}\n"
	if [[ "${k27file}" == "yes" ]]; then	
		for sample in ${tf_sample_list[@]} H3K27ac
		do
			printf "${sample}\n" > combined/peaks/temp_col_${analysisname}_${sample}.txt
			awk -v OFS="\t" -v s=${sample} '{if ($0 ~ s) print "1"; else print "0"}' combined/peaks/TF_peaks_${analysisname}.bed >> combined/peaks/temp_col_${analysisname}_${sample}.txt
		done
	else
		for sample in ${tf_sample_list[@]}
		do
			printf "${sample}\n" > combined/peaks/temp_col_${analysisname}_${sample}.txt
			awk -v OFS="\t" -v s=${sample} '{if ($0 ~ s) print "1"; else print "0"}' combined/peaks/TF_peaks_${analysisname}.bed >> combined/peaks/temp_col_${analysisname}_${sample}.txt
		done
	fi
	#### To group peaks based on their distance (gene body (x=0), promoter (0<x<2kb upstream), terminator (0<x<2kb downstream), distal)
	awk -v OFS="\t" 'BEGIN {printf "PeakID\tDistance\tGroup\n"} {if ($5<-2000) {d="Distal_downstream"; a=-$5} else if ($5<0) {d="Terminator"; a=-$5} else if ($5==0) {d="Gene_body"; a=$5} else if ($5>2000) {d="Distal_upstream"; a=$5} else {d="Promoter"; a=$5} print $4,a,d}' combined/peaks/TF_peaks_${analysisname}.bed > combined/peaks/temp_col_${analysisname}_AAA.txt
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
		# H3K4me1|H3K27me1|H3K27me2|H3K27me3|H3K9me1|H3K9me2|H3K9me3) 	peaktype="broad"
					# l=500
					# g=400;;
		# H3K27ac|H3K4me3) 	peaktype="narrow"
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

if [[ "${total}" == "NO" ]]; then
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

if [[ "${total}" != "TEST" ]]; then
	printf "\nDoing complete analysis for ${analysisname} with deeptools version:\n"
	deeptools --version
else
	printf "\nDoing testing analysis for ${analysisname}\n"
fi

uniq_chip_mark_list=($(printf "%s\n" "${chip_mark_list[@]}" | sort -u))

if [[ "${total}" != "TEST" ]]; then
	#### Splitting the region file by strand
	awk -v OFS="\t" '$6=="+"' ${regionfile} > combined/matrix/temp_regions_${regionname}_plus.bed
	awk -v OFS="\t" '$6=="-"' ${regionfile} > combined/matrix/temp_regions_${regionname}_minus.bed
	regionlabel=$(wc -l ${regionfile} | awk -v n=${regionname} '{print n"("$1")"}')

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
		computeMatrix scale-regions -q --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${regionname}_${strand}.bed -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/temp_all_genes_regions_${analysisname}_${strand}.gz
		printf "\nComputing reference-point on TSS ${strand} strand matrix for ${analysisname}\n"
		computeMatrix reference-point --referencePoint "TSS" -q --missingDataAsZero --skipZeros -R combined/matrix/temp_regions_${regionname}_${strand}.bed -S ${bw_list} -bs 50 -b 2000 -a 8000 -p ${threads} -o combined/matrix/temp_all_genes_tss_${analysisname}_${strand}.gz
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
		computeMatrixOperations rbind -m combined/matrix/temp_all_genes_${matrix}_${analysisname}_plus.gz combined/matrix/temp_all_genes_${matrix}_${analysisname}_minus.gz -o combined/matrix/all_genes_${matrix}_${analysisname}.gz
		printf "\nGetting scales for ${matrix} matrix of ${analysisname}\n"
		computeMatrixOperations dataRange -m combined/matrix/all_genes_${matrix}_${analysisname}.gz > combined/matrix/temp_values_${matrix}_${analysisname}.txt
		plotProfile -m combined/matrix/all_genes_${matrix}_${analysisname}.gz -out combined/plots/temp_${matrix}_${analysisname}_profile.pdf --samplesLabel ${all_labels[@]} --averageType mean --outFileNameData combined/matrix/temp_values_profile_${matrix}_${analysisname}.txt
		rm -f combined/plots/temp_${matrix}_${analysisname}_profile.pdf
		mins=()
		maxs=()
		ymins=()
		ymaxs=()
		for mark in ${all_samples[@]}
		do
			mini=$(grep "${mark}" combined/matrix/temp_values_${matrix}_${analysisname}.txt | awk 'BEGIN {m=999999} {a=$5; if (a<m) m=a;} END {print m}')
			maxi=$(grep "${mark}" combined/matrix/temp_values_${matrix}_${analysisname}.txt | awk 'BEGIN {m=-999999} {a=$6; if (a>m) m=a;} END {print m}')
			num=$(grep "${mark}" combined/matrix/temp_values_${matrix}_${analysisname}.txt | wc -l)
			test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
			if [[ ${test} == "yes" ]]; then
				mini=("0")
				maxi=("0.005")
			fi
			for i in $(seq 1 ${num})
			do
				mins+=("${mini}")
				maxs+=("${maxi}")
			done		
			ymini=$(grep "${mark}" combined/matrix/temp_values_profile_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
			ymaxi=$(grep "${mark}" combined/matrix/temp_values_profile_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
			num=$(grep "${mark}" combined/matrix/temp_values_profile_${matrix}_${analysisname}.txt | wc -l)
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
			mini=$(grep ${sample} combined/matrix/temp_values_${matrix}_${analysisname}.txt | awk '{print $5}')
			maxi=$(grep ${sample} combined/matrix/temp_values_${matrix}_${analysisname}.txt | awk '{print $6}')
			test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
			if [[ ${test} == "yes" ]]; then
				mins2+=("0")
				maxs2+=("0.005")
			else
				mins2+=("${mini}")
				maxs2+=("${maxi}")
			fi
		done
		ymins2=()
		ymaxs2=()
		for sample in ${all_labels[@]}
		do
			ymini=$(grep ${sample} combined/matrix/temp_values_profile_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
			ymaxi=$(grep ${sample} combined/matrix/temp_values_profile_${matrix}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
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
		plotHeatmap -m combined/matrix/all_genes_${matrix}_${analysisname}.gz -out combined/plots/all_genes_${analysisname}_heatmap_${matrix}.pdf --sortRegions descend --sortUsing mean --samplesLabel ${all_labels[@]} --regionsLabel ${regionlabel} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
		printf "\nPlotting heatmap for ${matrix} matrix of ${analysisname} scaling by sample\n"
		plotHeatmap -m combined/matrix/all_genes_${matrix}_${analysisname}.gz -out combined/plots/all_genes_${analysisname}_heatmap_${matrix}_v2.pdf --sortRegions descend --sortUsing mean --samplesLabel ${all_labels[@]} --regionsLabel ${regionlabel} --colorMap 'seismic' --zMin ${mins2[@]} --zMax ${maxs2[@]} --yMin ${ymins2[@]} --yMax ${ymaxs2[@]} --interpolationMethod 'bilinear'
	done
	rm -f combined/matrix/temp*${analysisname}*
fi

#########################################################################################
####################################### PART9 ###########################################
############################## Making heatmaps on DEGs  #################################
#########################################################################################

#### To make heatmaps and profiles with deeptools on the DEG if they were called

if [ ${#rnaseq_name_list[@]} -ge 2 ] && [[ "${total}" != "TEST" ]]; then
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
				computeMatrixOperations dataRange -m combined/matrix/${analysisname}_only_${namei}_DEG.gz > combined/matrix/temp_values_${analysisname}_only_${namei}_DEG.txt
				mins=()
				maxs=()
				for sample in ${sorted_labels[@]}
				do
					mini=$(grep ${sample} combined/matrix/temp_values_${analysisname}_only_${namei}_DEG.txt | awk '{print $5}')
					mins+=("${mini}")
					maxi=$(grep ${sample} combined/matrix/temp_values_${analysisname}_only_${namei}_DEG.txt | awk '{print $6}')
					maxs+=("${maxi}")
				done
				plotProfile -m combined/matrix/${analysisname}_only_${namei}_DEG.gz -out combined/plots/regions_${analysisname}_profile_only_${namei}_DEG.pdf --plotType 'lines' --averageType 'mean' --samplesLabel ${sorted_labels[@]} --regionsLabel "${namei}_UP" "${namei}_DOWN" --outFileNameData combined/matrix/temp2_values_${analysisname}_profile_only_${namei}_DEG.txt
				ymins=()
				ymaxs=()
				for sample in ${sorted_labels[@]}
				do
			 		ymini=$(grep ${sample} combined/matrix/temp2_values_${analysisname}_profile_only_${namei}_DEG.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
					ymaxi=$(grep ${sample} combined/matrix/temp2_values_${analysisname}_profile_only_${namei}_DEG.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
					ymins+=("${ymini}")
					ymaxs+=("${ymaxi}")
				done
				printf "\nPlotting heatmap for ${namei} specific DEG from ${analysisname}\n"
				plotHeatmap -m combined/matrix/${analysisname}_only_${namei}_DEG.gz -out combined/plots/${analysisname}_heatmap_only_${namei}_DEG.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sorted_labels[@]} --regionsLabel "${namei}_UP(${nbup})" "${namei}_DOWN(${nbdown})" --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --colorMap 'seismic' --interpolationMethod 'bilinear'
				plotProfile -m combined/matrix/${analysisname}_only_${namei}_DEG.gz -out combined/plots/regions_${analysisname}_profile_only_${namei}_DEG.pdf --plotType 'lines' --averageType 'mean' --samplesLabel ${sorted_labels[@]} --regionsLabel "${namei}_UP(${nbup})" "${namei}_DOWN(${nbdown})" --yMin ${ymins[@]} --yMax ${ymaxs[@]}
			done
		fi
	fi
	rm -f combined/matrix/*${analysisname}*only*
	rm -f combined/matrix/temp*${analysisname}*
fi


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
		awk -v OFS="\t" -v d="${cols}" -v t=${reps} 'BEGIN {split(d, a, " ")} NR > 1 {b=0; for (i in a) b+=$(a[i]); c=b/t; print $1,c}' combined/DEG/counts_${analysisname}.txt > combined/DEG/temp_counts_${analysisname}_${tissue}.txt
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
		if [[ "${total}" != "TEST" ]]; then
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
				computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/temp_split_expression_regions_${analysisname}_${strand}.gz --quiet
				printf "\nComputing reference-point on TSS $strand strand matrix for split expression in ${tissue} in $analysisname\n"
				computeMatrix reference-point --referencePoint "TSS" --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 8000 -p ${threads} -o combined/matrix/temp_split_expression_tss_${analysisname}_${strand}.gz --quiet
			done

			### Merging stranded matrix, extracting scales and plotting heatmaps
			for matrix in regions tss
			do
				tissue_labels="${tissue_labels_chip[*]} ${tissue_labels_rna[*]} ${tissue_labels_rampage[*]} ${tissue_labels_shrna[*]}"
				printf "\nMerging stranded matrices aligned by ${matrix} for split expression in ${tissue} in ${analysisname}\n"
				computeMatrixOperations rbind -m combined/matrix/temp_split_expression_${matrix}_${analysisname}_plus.gz combined/matrix/temp_split_expression_${matrix}_${analysisname}_minus.gz -o combined/matrix/temp_split_expression_${matrix}_${analysisname}.gz
				printf "\nGetting scales for ${matrix} matrix for split expression in ${tissue} in ${analysisname}\n"
				computeMatrixOperations dataRange -m combined/matrix/temp_split_expression_${matrix}_${analysisname}.gz > combined/matrix/temp_values_split_expression_${matrix}_${analysisname}.txt
				mins=()
				maxs=()
				for sample in ${tissue_labels[@]}
				do
					mini=$(grep ${sample} combined/matrix/temp_values_split_expression_${matrix}_${analysisname}.txt | awk '{print $5}')
					maxi=$(grep ${sample} combined/matrix/temp_values_split_expression_${matrix}_${analysisname}.txt | awk '{print $6}')
					test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
					if [[ ${test} == "yes" ]]; then
						mins+=("0")
						maxs+=("0.005")
					else
						mins+=("${mini}")
						maxs+=("${maxi}")
					fi
				done
				computeMatrixOperations sort -m combined/matrix/temp_split_expression_${matrix}_${analysisname}.gz -R ${sorted_regions[@]} -o combined/matrix/split_expression_${matrix}_${analysisname}.gz
				plotProfile -m combined/matrix/split_expression_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_profile_${matrix}.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --outFileNameData combined/matrix/temp2_values_split_expression_${matrix}_${tissue}_${analysisname}.txt
				ymins=()
				ymaxs=()
				for sample in ${tissue_labels[@]}
				do
		 			ymini=$(grep ${sample} combined/matrix/temp2_values_split_expression_${matrix}_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
					ymaxi=$(grep ${sample} combined/matrix/temp2_values_split_expression_${matrix}_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
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
				plotHeatmap -m combined/matrix/split_expression_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_heatmap_${matrix}.pdf --sortRegions keep --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
				if [[ ${matrix} == "tss" ]]; then
					printf "\nPlotting heatmap for ${matrix} matrix for split expression in ${tissue} in ${analysisname} scaling by sample by sorting by length\n"
					plotHeatmap -m combined/matrix/split_expression_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_heatmap_${matrix}_v2.pdf --sortRegions descend --sortUsing region_length --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
				fi
				printf "\nPlotting profile for ${matrix} matrix for split expression in ${tissue} in ${analysisname} scaling by sample\n"
				plotProfile -m combined/matrix/split_expression_${matrix}_${analysisname}.gz -out combined/plots/split_expression_${tissue}_${analysisname}_profile_${matrix}.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --yMin ${ymins[@]} --yMax ${ymaxs[@]}
			done
			rm -f combined/matrix/temp*${analysisname}*
		fi
		rm -f combined/DEG/temp*${analysisname}*
	fi
done

#########################################################################################
####################################### PART11 ##########################################
############ Making heatmaps on distal H3K27ac peaks split by ChIP enrichment  ##########
#########################################################################################

#### To make heatmap and profile with deeptools for each tissue based on grouped H3K27ac levels at distal elements (>2kb)

uniq_chip_tissue_list=($(printf "%s\n" "${chip_tissue_list[@]}" | sort -u))

h3k27actissues=()
createfile="no"
enhancerfile="no"
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
		printf "\nMaking heatmaps of distal enhancers (H3K27ac peak >2kb from TSS) in ${tissue}\n"
		printf "\nGetting bed file of distal enhancers for ${tissue}\n"
		enhancerfile="yes"
		bedtools sort -g ${ref_dir}/chrom.sizes -i ChIP/peaks/best_peaks_${line}_${tissue}_H3K27ac.bed | awk '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/)'> combined/peaks/temp_${analysisname}_${tissue}.bed
		if [ -s combined/DEG/sorted_expression_${analysisname}_${tissue}.bed ]; then
			bedtools sort -g ${ref_dir}/chrom.sizes -i combined/DEG/sorted_expression_${analysisname}_${tissue}.bed > combined/peaks/temp_${analysisname}_${tissue}_expression.bed
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>= 2000 && $16=="+") print $1,$2+$10,$12,$4,$5,$16,$14,$15; else if ($17<= -2000 && $16=="-") print $1,$13,$2+$10,$4,$5,$16,$14,$15}' | sort -k5,5nr > combined/peaks/distal_${analysisname}_${tissue}.bed
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>= 2000 && $16=="+") print $1,$2,$3,$4,$5,$16,$14,$15; else if ($17<= -2000 && $16=="-") print $1,$2,$3,$4,$5,$16,$14,$15}' > combined/peaks/enhancers_distal_upstream_${analysisname}_${tissue}.txt
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>= 2000 && $16=="-") print $1,$2,$3,$4,$5,$16,$14,$15; else if ($17<= -2000 && $16=="+") print $1,$2,$3,$4,$5,$16,$14,$15}' > combined/peaks/enhancers_distal_downstream_${analysisname}_${tissue}.txt
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>0 && $17<2000 && $16=="+") print $1,$2,$3,$4,$5,$16,$14,$15; else if ($17>-2000 && $17<0 && $16=="-") print $1,$2,$3,$4,$5,$16,$14,$15}' > combined/peaks/enhancers_promoter_${analysisname}_${tissue}.txt
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>0 && $17<2000 && $16=="-") print $1,$2,$3,$4,$5,$16,$14,$15; else if ($17>-2000 && $17<0 && $16=="+") print $1,$2,$3,$4,$5,$16,$14,$15}' > combined/peaks/enhancers_terminator_${analysisname}_${tissue}.txt
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17==0) print $1,$2,$3,$4,$5,$16,$14,$15}' > combined/peaks/enhancers_genic_${analysisname}_${tissue}.txt
		else
			bedtools sort -g ${ref_dir}/chrom.sizes -i ${regionfile} > combined/peaks/temp_${analysisname}_${tissue}_no_expression.bed
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_no_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>= 2000 && $16=="+") print $1,$2+$10,$12,$4,$5,$16,$14; else if ($17<= -2000 && $16=="-") print $1,$13,$2+$10,$4,$5,$16,$14}' | sort -k5,5nr > combined/peaks/distal_${analysisname}_${tissue}.bed
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_no_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>= 2000 && $16=="+") print $1,$2,$3,$4,$5,$16,$14; else if ($17<= -2000 && $16=="-") print $1,$2,$3,$4,$5,$16,$14}' > combined/peaks/enhancers_distal_upstream_${analysisname}_${tissue}.txt
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_no_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>= 2000 && $16=="-") print $1,$2,$3,$4,$5,$16,$14; else if ($17<= -2000 && $16=="+") print $1,$2,$3,$4,$5,$16,$14}' > combined/peaks/enhancers_distal_downstream_${analysisname}_${tissue}.txt
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_no_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>0 && $17<2000 && $16=="+") print $1,$2,$3,$4,$5,$16,$14; else if ($17>-2000 && $17<0 && $16=="-") print $1,$2,$3,$4,$5,$16,$14}' > combined/peaks/enhancers_promoter_${analysisname}_${tissue}.txt
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_no_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17>0 && $17<2000 && $16=="-") print $1,$2,$3,$4,$5,$16,$14; else if ($17>-2000 && $17<0 && $16=="+") print $1,$2,$3,$4,$5,$16,$14}' > combined/peaks/enhancers_terminator_${analysisname}_${tissue}.txt
			bedtools closest -a combined/peaks/temp_${analysisname}_${tissue}.bed -b combined/peaks/temp_${analysisname}_${tissue}_no_expression.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{if ($17==0) print $1,$2,$3,$4,$5,$16,$14}' > combined/peaks/enhancers_genic_${analysisname}_${tissue}.txt
		fi
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
			if [ -s combined/DEG/sorted_expression_${analysisname}_${tissue}.bed ]; then
				awk -v OFS="\t" -v g=${name} -v t=${tissue} '{print t,$4,$5,$7,$8,g }' combined/peaks/temp_distal_${analysisname}_${tissue}_group${i}.bed >> combined/peaks/temp2_distal_peaks_${analysisname}_${tissue}.txt
				createfile="yes"
			fi
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
			printf "\nComputing scale-regions ${strand} strand matrix for ${tissue}\n"
			computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${sorted_regions[@]} -S ${bw_list} -bs 50 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/temp_distal_${tissue}_${analysisname}_${strand}.gz --quiet
		done
		printf "\nMerging stranded matrices for distal peaks in ${tissue}\n"
		computeMatrixOperations rbind -m combined/matrix/temp_distal_${tissue}_${analysisname}_plus.gz combined/matrix/temp_distal_${tissue}_${analysisname}_minus.gz -o combined/matrix/temp_distal_${tissue}_${analysisname}.gz
		printf "\nGetting scales for distal peaks in ${tissue}\n"
		computeMatrixOperations dataRange -m combined/matrix/temp_distal_${tissue}_${analysisname}.gz > combined/matrix/temp_values_distal_${tissue}_${analysisname}.txt
		mins=()
		maxs=()
		for sample in ${tissue_labels[@]}
		do
			mini=$(grep ${sample} combined/matrix/temp_values_distal_${tissue}_${analysisname}.txt | awk '{print $5}')
			maxi=$(grep ${sample} combined/matrix/temp_values_distal_${tissue}_${analysisname}.txt | awk '{print $6}')
			test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
			if [[ ${test} == "yes" ]]; then
				mins+=("0")
				maxs+=("0.005")
			else
				mins+=("${mini}")
				maxs+=("${maxi}")
			fi
		done
		computeMatrixOperations sort -m combined/matrix/temp_distal_${tissue}_${analysisname}.gz -R ${sorted_regions[@]} -o combined/matrix/distal_${tissue}_${analysisname}.gz
		plotProfile -m combined/matrix/distal_${tissue}_${analysisname}.gz -out combined/plots/distal_${tissue}_${analysisname}_profile_mean.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --outFileNameData combined/matrix/temp2_values_distal_${tissue}_${analysisname}.txt
		ymins=()
		ymaxs=()
		for sample in ${tissue_labels[@]}
		do
		 	ymini=$(grep ${sample} combined/matrix/temp2_values_distal_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
			ymaxi=$(grep ${sample} combined/matrix/temp2_values_distal_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
			test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
			if [[ ${test} == "yes" ]]; then
				ymins+=("0")
				ymaxs+=("0.01")
			else
				ymins+=("${ymini}")
				ymaxs+=("${ymaxi}")
			fi
		done
		printf "\nPlotting heatmap for distal peaks in ${tissue} in ${analysisname} scaling by sample\n"
		plotHeatmap -m combined/matrix/distal_${tissue}_${analysisname}.gz -out combined/plots/distal_${tissue}_${analysisname}_heatmap_mean.pdf --sortRegions keep --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --startLabel "enhancer" --endLabel "TSS"
		printf "\nPlotting mean profile for distal peaks in ${tissue} in ${analysisname} scaling by sample\n"
		plotProfile -m combined/matrix/distal_${tissue}_${analysisname}.gz -out combined/plots/distal_${tissue}_${analysisname}_profile_mean.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType mean --yMin ${ymins[@]} --yMax ${ymaxs[@]} --startLabel "enhancer" --endLabel "TSS"
		
		plotProfile -m combined/matrix/distal_${tissue}_${analysisname}.gz -out combined/plots/distal_${tissue}_${analysisname}_profile_median.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType median --outFileNameData combined/matrix/temp3_values_distal_${tissue}_${analysisname}.txt
		ymins=()
		ymaxs=()
		for sample in ${tissue_labels[@]}
		do
		 	ymini=$(grep ${sample} combined/matrix/temp3_values_distal_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
			ymaxi=$(grep ${sample} combined/matrix/temp3_values_distal_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
			test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
			if [[ ${test} == "yes" ]]; then
				ymins+=("0")
				ymaxs+=("0.01")
			else
				ymins+=("${ymini}")
				ymaxs+=("${ymaxi}")
			fi
		done
		printf "\nPlotting median profile for distal peaks in ${tissue} in ${analysisname} scaling by sample\n"
		plotProfile -m combined/matrix/distal_${tissue}_${analysisname}.gz -out combined/plots/distal_${tissue}_${analysisname}_profile_median.pdf --samplesLabel ${tissue_labels[@]} --regionsLabel ${regions_labels[@]} --averageType median --yMin ${ymins[@]} --yMax ${ymaxs[@]} --startLabel "enhancer" --endLabel "TSS"		
	else
		printf "\nTissue ${tissue} will not be processed\n"
	fi
done
if [[ "${createfile}" == "yes" ]]; then
	printf "Tissue\tPeak_ID\tPeakQuality\tGID\tRPKM\tGroup\n" > combined/peaks/all_grouped_distal_peaks_${analysisname}.txt
	cat combined/peaks/temp2_distal_peaks_${analysisname}_*.txt > combined/peaks/all_grouped_distal_peaks_${analysisname}.txt
fi

if [[ "${enhancerfile}" == "yes" ]]; then
	printf "Line\tTissue\tEnhancer\tCount\n" > combined/peaks/summary_enhancers_${analysisname}.txt
	for tissue in ${h3k27actissues[@]}
	do
		for type in genic promoter terminator distal_upstream distal_downstream 
		do
			wc -l combined/peaks/enhancers_${type}_${analysisname}_${tissue}.txt | awk -v OFS="\t" -v l=$line -v t=$tissue -v d=$type '{print l,t,d,$1}' >> combined/peaks/summary_enhancers_${analysisname}.txt
		done
	done
fi
rm -f combined/matrix/temp*${analysisname}*
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
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_RNAseq_plus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=1000*($5/($3-$2)); print a}' > combined/peaks/col_B_${analysisname}_${line}_${tissue}.txt
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
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_RNAseq_minus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=1000*($5/($3-$2)); print a}' > combined/peaks/col_C_${analysisname}_${line}_${tissue}.txt
			rm -f combined/peaks/temp_RNAseq_minus_${analysisname}_${line}.bedGraph
			header="${header}\tRNAseq_minus"
		fi
	done
	for bw in ${rampage_bw_list_plus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			printf "\nGetting RAMPAGE plus strand coverage on H3K27ac peaks for ${tissue}\n"
			bigWigToBedGraph ${bw} combined/peaks/temp_RAMPAGE_plus_${analysisname}_${line}.bedGraph
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_RAMPAGE_plus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=1000*($5/($3-$2)); print a}' > combined/peaks/col_D_${analysisname}_${line}_${tissue}.txt
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
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_RAMPAGE_minus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=1000*($5/($3-$2)); print a}' > combined/peaks/col_E_${analysisname}_${line}_${tissue}.txt
			rm -f combined/peaks/temp_RAMPAGE_minus_${analysisname}_${line}.bedGraph
			header="${header}\tRAMPAGE_minus"
		fi
	done
	for bw in ${shrna_bw_list_plus[*]}
	do					
		if [[ ${bw} =~ ${tissue} ]]; then
			printf "\nGetting shRNA plus strand coverage on H3K27ac peaks for ${tissue}\n"
			bigWigToBedGraph ${bw} combined/peaks/temp_shRNA_plus_${analysisname}_${line}.bedGraph
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_shRNA_plus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=1000*($5/($3-$2)); print a}' > combined/peaks/col_F_${analysisname}_${line}_${tissue}.txt
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
			bedtools intersect -a combined/peaks/col_A_${analysisname}_${line}_${tissue}.txt -b combined/peaks/temp_shRNA_minus_${analysisname}_${line}.bedGraph -wao | awk -v OFS="\t" '{if ($8 == ".") $5=0; else $5=$8*$9; print $1,$2,$3,$4,$5}' | bedtools merge -i stdin -o distinct,sum -c 4,5 | awk -v OFS="\t" '{a=1000*($5/($3-$2)); print a}' > combined/peaks/col_G_${analysisname}_${line}_${tissue}.txt
			rm -f combined/peaks/temp_shRNA_minus_${analysisname}_${line}.bedGraph
			header="${header}\tshRNA_minus"
		fi
	done	
	printf "${header}\n" > combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt
	paste combined/peaks/col_*_${analysisname}_${line}_${tissue}.txt >> combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt
	plot="No"
	if [[ ${rnaseq} == 1 ]] && [[ ${rampage} == 1 ]] && [[ ${shrna} == 1 ]]; then
		included_samples="RNAseq_and_RAMPAGE_and_shRNA"
		plot="Yes"
	elif [[ ${rnaseq} == 1 ]] && [[ ${rampage} == 1 ]]; then
		included_samples="RNAseq_and_RAMPAGE"
		plot="Yes"
	elif [[ ${rnaseq} == 1 ]] && [[ ${shrna} == 1 ]]; then
		included_samples="RNAseq_and_shRNA"
		plot="Yes"
	elif [[ ${rampage} == 1 ]] && [[ ${shrna} == 1 ]]; then
		included_samples="RAMPAGE_and_shRNA"
		plot="Yes"
	elif [[ ${rnaseq} == 1 ]]; then
		included_samples="RNAseq"
		plot="Yes"
	elif [[ ${rampage} == 1 ]]; then
		included_samples="RAMPAGE"
		plot="Yes"
	elif [[ ${shRNA} == 1 ]]; then
		included_samples="shRNA"
		plot="Yes"
	fi
	if [[ ${plot} == "Yes" ]]; then
		#### To plot correlation between RNA expression datasets at distal peaks
		printf "\nCreating scatter plot for ${analysisname} ${tissue} with R version:\n"
		R --version
		Rscript --vanilla ${mc_dir}/MaizeCode_R_scatter_distal_peaks.r ${analysisname} ${tissue} ${line} ${included_samples} combined/peaks/all_grouped_distal_peaks_${analysisname}.txt combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt
	fi
done

############################################################################################
######################################### PART13 ###########################################
############# Gathering all available information on each type of enhancers ################
############################################################################################

rnaseq="No"
rampage="No"
shrna="No"
deg="No"
tf="No"
for tissue in ${h3k27actissues[@]}
do
	if [ $(grep "RNAseq" combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt | wc -l) -gt 0 ]; then
		rnaseq="Yes"
	fi
	if [ $(grep "RAMPAGE" combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt | wc -l) -gt 0 ]; then
		rampage="Yes"
	fi
	if [ $(grep "shRNA" combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt | wc -l) -gt 0 ]; then
		shrna="Yes"
	fi
	if [ -s combined/peaks/all_${line}_${tissue}_${analysisname}_DEG_GID.txt ]; then
		deg="Yes"
	fi
	if [ -s combined/peaks/TF_peaks_${analysisname}.bed ]; then
		tf="Yes"
	fi
	for type in distal_upstream distal_downstream promoter terminator genic
	do
		printf "Gathering all available RNA (${rnaseq}), RAMPAGE (${rampage}), shRNA (${shrna}), DEG (${deg}) and TF (${tf}) data for ${type} enhancers of ${line} ${tissue}\n"
		if [ -s combined/peaks/temp_complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt ]; then
			rm -f combined/peaks/temp_complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
		fi
		while read chr start stop peakID quality strand GID expression
		do
			colnb=7
			header="Chr\tStart\tStop\tPeakID\tQuality\tstrand\tGID"
			rowi="${chr}\t${start}\t${stop}\t${peakID}\t${quality}\t${strand}\t${GID}"
			if [[ ${rnaseq} == "Yes" ]]; then	
				RNAseq_plus=$(awk -v p=${peakID} '$4 == p {print $5}' combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt)
				RNAseq_minus=$(awk -v p=${peakID} '$4 == p {print $6}' combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt)
				rowi="${rowi}\t${expression}\t${RNAseq_plus}\t${RNAseq_minus}"
				header="${header}\texpression\tRNAseq_plus\tRNAseq_minus"
				colnb=$((colnb+3))
			fi
			if [[ ${rampage} == "Yes" ]]; then
				RAMPAGE_plus=$(awk -v p=${peakID} '$4 == p {print $7}' combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt)
				RAMPAGE_minus=$(awk -v p=${peakID} '$4 == p {print $8}' combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt)
				header="${header}\tRAMPAGE_plus\tRAMPAGE_minus"
				rowi="${rowi}\t${RAMPAGE_plus}\t${RAMPAGE_minus}"
				colnb=$((colnb+2))
			fi
			if [[ ${shrna} == "Yes" ]]; then
				shRNA_plus=$(awk -v p=${peakID} '$4 == p {print $9}' combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt)
				shRNA_minus=$(awk -v p=${peakID} '$4 == p {print $10}' combined/peaks/H3K27ac_peaks_expression_${line}_${tissue}_${analysisname}.txt)
				header="${header}\tshRNA_plus\tshRNA_minus"
				rowi="${rowi}\t${shRNA_plus}\t${shRNA_minus}"
				colnb=$((colnb+2))
			fi
			if [[ ${deg} == "Yes" ]]; then
				if grep -q "$GID" combined/DEG/only_${line}_${tissue}_DEG_UP_${analysisname}.bed
				then
					DEG="unique_UP"
				elif grep -q "$GID" combined/DEG/only_${line}_${tissue}_DEG_DOWN_${analysisname}.bed
				then
					DEG="unique_DOWN"
				elif grep -q "$GID" combined/peaks/all_${line}_${tissue}_${analysisname}_DEG_GID.txt
				then
					DEG="DEG"
				else
					DEG="Not_DEG"
				fi
				header="${header}\tDEG"
				rowi="${rowi}\t${DEG}"
				colnb=$((colnb+1))
			fi
			printf "${rowi}\n" >> combined/peaks/temp_complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
		done < combined/peaks/enhancers_${type}_${analysisname}_${tissue}.txt
		if [[ ${tf} == "Yes" ]]; then
			printf "${header}\tTFs\n" > combined/peaks/complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
			array="4"
			limit=$((colnb+1))
			for (( i=5; i<=${limit}; i++ ))
			do
				array="${array},${i}"
			done
			bedtools sort -g ${ref_dir}/chrom.sizes -i combined/peaks/temp_complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt > combined/peaks/temp2_complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
			bedtools sort -g ${ref_dir}/chrom.sizes -i combined/peaks/TF_peaks_${analysisname}.bed > combined/peaks/temp3_TFs_${analysisname}.txt
			tfcol=$((colnb+7))
			bedtools intersect -wao -a combined/peaks/temp2_complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt -b combined/peaks/temp3_TFs_${analysisname}.txt | awk -v OFS="\t" -v c=${tfcol} -v a=${colnb} '{if ($c == ".") t="None"; else t=$c ; for (i=1; i<=a; i++) {printf $i"\t"}; print t}' > combined/peaks/temp4_complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
			bedtools merge -o distinct -c ${array} -i combined/peaks/temp4_complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt >> combined/peaks/complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
		else
			printf "${header}\n" > combined/peaks/complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
			bedtools sort -g ${ref_dir}/chrom.sizes -i combined/peaks/temp_complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt >> combined/peaks/complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
		fi
		rm -f combined/peaks/temp*${analysisname}*
	done
done

############################################################################################
######################################### PART14 ###########################################
####################### Making heatmaps on all types of enhancers ##########################
############################################################################################

awk -v OFS="\t" '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/) {print $1,$2,$3,"1"}' ${regionfile} | bedtools sort -g ${ref_dir}/chrom.sizes > ChIP/tracks/temp_${regionname}.bg
bedtools merge -i ChIP/tracks/temp_${regionname}.bg -o max -c 4 | LC_COLLATE=C sort -k1,1 -k2,2n > ChIP/tracks/temp2_${regionname}.bg
bedGraphToBigWig ChIP/tracks/temp2_${regionname}.bg ${ref_dir}/chrom.sizes ChIP/tracks/${regionname}.bw
rm -f ChIP/tracks/temp*.bg

### Need to add the TE file for all other inbreds (based on reference name)
### Maybe changing the way to use the references altogether, like adding a reference folder in the git hub and have them add theirs there

tefilebw=""
tefilebed=""
if [ -s /grid/martienssen/data_norepl/dropbox/maizecode/TEs/${ref}_TEs.gff3.gz ]; then
	if [ ! -d combined/TSS ]; then
		mkdir combined/TSS
	fi	
	if [ ! -s combined/TSS/${ref}_all_tes.bed ]; then
		zcat /grid/martienssen/data_norepl/dropbox/maizecode/TEs/${ref}_TEs.gff3.gz | awk -v OFS="\t" '$1 !~ /^#/ {print $1,$4-1,$5,$3,".",$7}' | bedtools sort -g ${ref_dir}/chrom.sizes > combined/TSS/${ref}_all_tes.bed
	fi
	awk -v OFS="\t" '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/) {print $1,$2,$3,"1"}' combined/TSS/${ref}_all_tes.bed | bedtools sort -g ${ref_dir}/chrom.sizes > ChIP/tracks/temp_${ref}_all_tes.bg
	bedtools merge -i ChIP/tracks/temp_${ref}_all_tes.bg -o max -c 4 | LC_COLLATE=C sort -k1,1 -k2,2n > ChIP/tracks/temp2_${ref}_all_tes.bg
	bedGraphToBigWig ChIP/tracks/temp2_${ref}_all_tes.bg ${ref_dir}/chrom.sizes ChIP/tracks/${ref}_all_tes.bw
	tefilebw="ChIP/tracks/${ref}_all_tes.bw"
	tefilebed="combined/TSS/${ref}_all_tes.bed"
	rm -f ChIP/tracks/temp*.bg
	awk '{print $4}' combined/TSS/${ref}_all_tes.bed | sort -u > combined/TSS/${ref}_TE_types.txt
	TEtypestring=""
	while read TEtype
	do
		TEtypestring="${TEtypestring},${TEtype}"
	done < combined/TSS/${ref}_TE_types.txt
fi

for tissue in ${h3k27actissues[@]}
do
	tissue_labels=()
	tissue_bw_plus=()
	tissue_bw_minus=()
	rnaseq=0
	rampage=0
	for sample in ${chip_sample_list[@]}
	do
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
			rnaseq=1
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
			rampage=1
		fi
	done
	if [[ ${rnaseq} == 1 ]] && [[ ${rampage} == 1 ]]; then
		for type in genic promoter terminator distal_upstream distal_downstream
		do
			awk -v OFS="\t" 'NR>1 {print $1,$2,$3,$4,$5,$6,$9+$10,$11+$12,$7}' combined/peaks/complete_enhancers_${type}_${line}_${tissue}_${analysisname}.txt | sort -k7,7nr -k8,8nr > combined/peaks/sorted_enhancers_${type}_${line}_${tissue}_${analysisname}.txt			
			awk '$6=="+"' combined/peaks/sorted_enhancers_${type}_${line}_${tissue}_${analysisname}.txt > combined/peaks/sorted_enhancers_${type}_${line}_${tissue}_${analysisname}_plus.txt
			awk '$6=="-"' combined/peaks/sorted_enhancers_${type}_${line}_${tissue}_${analysisname}.txt > combined/peaks/sorted_enhancers_${type}_${line}_${tissue}_${analysisname}_minus.txt
			nb=$(wc -l combined/peaks/sorted_enhancers_${type}_${line}_${tissue}_${analysisname}.txt | awk '{print $1}')
			regions_label="${type}(${nb})"
			
			for strand in plus minus
			do
				case "$strand" in
					plus)	bw="${tissue_bw_plus[*]} ChIP/tracks/${regionname}.bw ${tefilebw}"
						regions="combined/peaks/sorted_enhancers_${type}_${line}_${tissue}_${analysisname}_plus.txt";;
					minus)	bw="${tissue_bw_minus[*]} ChIP/tracks/${regionname}.bw ${tefilebw}"
						regions="combined/peaks/sorted_enhancers_${type}_${line}_${tissue}_${analysisname}_minus.txt";;
				esac
				printf "\nComputing matrix for ${line} ${tissue} ${type} ${strand} strand\n"
				computeMatrix scale-regions -q --missingDataAsZero --skipZeros -R ${regions} -S ${bw[@]} -bs 10 -b 3000 -a 3000 -m 1000 -p ${threads} -o combined/matrix/temp_regions_enhancers_${type}_${line}_${tissue}_${analysisname}_${strand}.gz
			done
			computeMatrixOperations rbind -m combined/matrix/temp_regions_enhancers_${type}_${line}_${tissue}_${analysisname}_plus.gz combined/matrix/temp_regions_enhancers_${type}_${line}_${tissue}_${analysisname}_minus.gz -o combined/matrix/temp_regions_enhancers_${type}_${line}_${tissue}_${analysisname}.gz
			if [[ ${tefilebw} != "" ]]; then
				label_list="${tissue_labels[*]} Genes TEs"
			else
				label_list="${tissue_labels[*]} Genes"
			fi
			printf "\nGetting scales for ${line} ${tissue} ${type}\n"
			plotProfile -m combined/matrix/temp_regions_enhancers_${type}_${line}_${tissue}_${analysisname}.gz -out combined/plots/enhancers_${type}_${line}_${tissue}_${analysisname}_temp_profile.pdf --samplesLabel ${label_list[@]} --averageType mean --outFileNameData combined/matrix/temp_values_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
			rm -f combined/plots/enhancers_${type}_${line}_${tissue}_${analysisname}_temp_profile.pdf
			ymins=()
			ymaxs=()
			for sample in ${label_list[@]}
			do
				ymini1=$(grep $sample combined/matrix/temp_values_enhancers_${type}_${line}_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.5; else a=m*0.5; print a}')
				ymini=$(printf "%.18f" "${ymini1}")
				ymaxi=$(grep $sample combined/matrix/temp_values_enhancers_${type}_${line}_${tissue}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m*1.5}')
				test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
				if [[ "${test}" == "yes" ]]; then
					ymins+=("0")
					ymaxs+=("0.01")
				else
					ymins+=("${ymini}")
					ymaxs+=("${ymaxi}")
				fi
			done
			computeMatrixOperations dataRange -m combined/matrix/temp_regions_enhancers_${type}_${line}_${tissue}_${analysisname}.gz > combined/matrix/temp2_values_enhancers_${type}_${line}_${tissue}_${analysisname}.txt
			mins=()
			maxs=()
			arr=${#tissue_labels[*]}
			for (( i=1; i<=${arr}; i++ ))
			do 
				mini=$(awk -v i=$i 'NR==(i+1) {print $5}' combined/matrix/temp2_values_enhancers_${type}_${line}_${tissue}_${analysisname}.txt)		
				maxi=$(awk -v i=$i 'NR==(i+1) {print $6}' combined/matrix/temp2_values_enhancers_${type}_${line}_${tissue}_${analysisname}.txt)
				test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
				if [[ "${test}" == "yes" ]]; then
					mins+=("0")
					maxs+=("0.005")
				else
					mins+=("${mini}")
					maxs+=("${maxi}")
				fi
			done
			if [[ ${tefilebw} != "" ]]; then
				mins+=("0" "0")
				maxs+=("1" "1")
			else
				mins+=("0")
				maxs+=("1")
			fi
			printf "\nSorting matrix for ${line} ${tissue} ${type}\n"
			computeMatrixOperations sort -m combined/matrix/temp_regions_enhancers_${type}_${line}_${tissue}_${analysisname}.gz -R combined/peaks/sorted_enhancers_${type}_${line}_${tissue}_${analysisname}.txt -o combined/matrix/sorted_regions_enhancers_${type}_${line}_${tissue}_${analysisname}.gz
			printf "\nPlotting heatmap for ${line} ${tissue} ${type}\n"
			plotHeatmap -m combined/matrix/sorted_regions_enhancers_${type}_${line}_${tissue}_${analysisname}.gz -out combined/plots/enhancers_${type}_${line}_${tissue}_${analysisname}_sortedbyRNA.pdf --sortRegions keep --samplesLabel ${label_list[*]} --regionsLabel ${regions_label} --colorMap 'seismic' --interpolationMethod 'bilinear' --yMin ${ymins[@]} --yMax ${ymaxs[@]} --zMin ${mins[@]} --zMax ${maxs[@]}
			plotHeatmap -m combined/matrix/sorted_regions_enhancers_${type}_${line}_${tissue}_${analysisname}.gz -out combined/plots/enhancers_${type}_${line}_${tissue}_${analysisname}_sortedbyAUTO.pdf --sortRegions descend --sortUsing mean --samplesLabel ${label_list[*]} --regionsLabel ${regions_label} --colorMap 'seismic' --interpolationMethod 'bilinear' --yMin ${ymins[@]} --yMax ${ymaxs[@]} --zMin ${mins[@]} --zMax ${maxs[@]}
		done
		rm -f combined/peaks/temp*${analysisname}*
		rm -f combined/matrix/temp*${analysisname}*
	fi
done

############################################################################################
########################################## PART14 ##########################################
############################## RAMPAGE TSS at genes and TEs ################################
############################################################################################

uniq_rampage_tissue_list=($(printf "%s\n" "${rampage_tissue_list[@]}" | sort -u))

if [[ ${#uniq_rampage_tissue_list[*]} -ge 2 ]] && [[ ${tefilebw} != "" ]]; then
	if [ -e combined/TSS/tmp_TSS_peaks_${analysisname}.bed ]; then
		rm -f combined/TSS/tmp_TSS_peaks_${analysisname}.bed
	fi
	for tissue in ${uniq_rampage_tissue_list[@]}
	do
		printf "\nMaking TSS peak file for ${tissue}\n"
		awk -v OFS="\t" '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {print $1,$2,$3,"Peak_"NR,"."}' RNA/TSS/idr_${line}_${tissue}_RAMPAGE.narrowPeak | bedtools sort -g ${ref_dir}/chrom.sizes > combined/TSS/${line}_${tissue}_RAMPAGE_peaks.bed
		awk -v OFS="\t" -v t=${tissue} '{print $1,$2,$3,t}' combined/TSS/${line}_${tissue}_RAMPAGE_peaks.bed | sort -k1,1 -k2,2n -u >> combined/TSS/tmp_TSS_peaks_${analysisname}.bed
		printf "\nGetting closest gene for ${tissue}\n"	
		bedtools closest -a combined/TSS/${line}_${tissue}_RAMPAGE_peaks.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/TSS/temp_TSS_${analysisname}.bed
		printf "\nGrouping based on distance for ${tissue}\n"
		awk -v OFS="\t" '{if ($5<-2000) {d="Intergenic"} else if ($5<0) {d="Terminator"} else if ($5==0) {d="Gene_body"} else if ($5>2000) {d="Intergenic"} else {d="Promoter"} print $0,d}' combined/TSS/temp_TSS_${analysisname}.bed > combined/TSS/temp2_TSS_${analysisname}.bed
		printf "\nIntersecting TE for ${tissue}\n"
		bedtools intersect -a combined/TSS/temp2_TSS_${analysisname}.bed -b combined/TSS/${ref}_all_tes.bed -loj | awk -v OFS="\t" -v t=${tissue} -v l=${line} '{if ($13==".") print l,t,l"_"t"_"$4,$9,"No",$9,$9; else if ($9 == "Intergenic") print l,t,l"_"t"_"$4,$9,$13,$13,$13; else print l,t,l"_"t"_"$4,$9,$13,$13,$13"_in_"$9}' > combined/TSS/TSS_in_genes_and_tes_${tissue}_${analysisname}.bed
		rm -f combined/TSS/temp*_TSS_${analysisname}.bed
	done
	printf "Line\tTissue\tPeak_ID\tGene\tTE\tLabel\tLabelcombined\n" > combined/TSS/Table_TSS_tissues_${analysisname}.txt
	cat combined/TSS/TSS_in_genes_and_tes_*_${analysisname}.bed >> combined/TSS/Table_TSS_tissues_${analysisname}.txt
	
	printf "\nPreparing merged TSS file for ${analysisname}\n"
	sort -k1,1 -k2,2n combined/TSS/tmp_TSS_peaks_${analysisname}.bed > combined/TSS/tmp2_TSS_peaks_${analysisname}.bed
	bedtools merge -i combined/TSS/tmp2_TSS_peaks_${analysisname}.bed -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{print $1,$2,$3,"Peak_"NR,$4}'> combined/TSS/tmp3_TSS_peaks_${analysisname}.bed
	printf "\nGetting closest gene for TSS in ${analysisname}\n"
	if [[ ${ref} == "B73_v4" ]]; then
		bedtools closest -a combined/TSS/tmp3_TSS_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/TSS/tmp4_TSS_peaks_${analysisname}.bed
	else
		bedtools closest -a combined/TSS/tmp3_TSS_peaks_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/TSS/tmp4_TSS_peaks_${analysisname}.bed
	fi
	printf "\nGrouping based on distance\n"
	awk -v OFS="\t" '{if ($5<-2000) {d="Intergenic"} else if ($5<0) {d="Terminator"} else if ($5==0) {d="Gene_body"} else if ($5>2000) {d="Intergenic"} else {d="Promoter"} print $0,d}' combined/TSS/tmp4_TSS_peaks_${analysisname}.bed > combined/TSS/tmp5_TSS_peaks_${analysisname}.bed
	printf "\nIntersecting with TEs\n"
	bedtools intersect -a combined/TSS/tmp5_TSS_peaks_${analysisname}.bed -b combined/TSS/${ref}_all_tes.bed -loj | awk -v OFS="\t" -v l=${line} 'BEGIN {printf "Line\tPeak_ID\tGene\tTE\tLabel\tLabelcombined\tGID\tTissues\tPeak_coordinates\n"} {if ($13==".") print l,$4,$9,"No",$9,$9,$8,$7,$1":"$2"-"$3; else if ($9 == "Intergenic") print l,$4,$9,$13,$13,$13,$8,$7,$1":"$2"-"$3; else print l,$4,$9,$13,$13,$13"_in_"$9,$8,$7,$1":"$2"-"$3}' > combined/TSS/all_TSS_in_genes_and_tes_${analysisname}.bed
	rm -f combined/TSS/tmp*_TSS_peaks_${analysisname}.bed
	#### To create a matrix of peak presence in each sample
	printf "\nCreating matrix file for ${analysisname}\n"
	for tissue in ${uniq_rampage_tissue_list[@]}
	do
		printf "${tissue}\n" > combined/TSS/temp_col_TSS_${analysisname}_${tissue}.txt
		awk -v OFS="\t" -v t=${tissue} 'NR>1 {if ($8 ~ t) print "1"; else print "0"}' combined/TSS/all_TSS_in_genes_and_tes_${analysisname}.bed >> combined/TSS/temp_col_TSS_${analysisname}_${tissue}.txt
	done
	awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7}' combined/TSS/all_TSS_in_genes_and_tes_${analysisname}.bed > combined/TSS/temp_col_TSS_${analysisname}_AAA.txt
	paste combined/TSS/temp_col_TSS_${analysisname}_*.txt | uniq > combined/TSS/matrix_upset_TSS_${analysisname}.txt
	rm -f combined/TSS/temp_col_TSS_${analysisname}_*.txt
	#### To make an Upset plot highlighting peaks in gene bodies
	printf "\nCreating distribution and Upset plot for TSS in ${analysisname} with R version:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_TSS_distribution_upset.r ${analysisname} ${TEtypestring} combined/TSS/Table_TSS_tissues_${analysisname}.txt combined/TSS/matrix_upset_TSS_${analysisname}.txt combined/TSS/all_TSS_in_genes_and_tes_${analysisname}.bed combined/DEG/genes_rpkm_${analysisname}.txt
fi

############################################################################################
########################################## PART15 ##########################################
############################ shRNA clusters at genes and TEs ###############################
############################################################################################

uniq_shrna_tissue_list=($(printf "%s\n" "${shrna_tissue_list[@]}" | sort -u))

if [[ ${#uniq_shrna_tissue_list[*]} -ge 2 ]] && [[ ${tefilebw} != "" ]]; then
	if [ ! -d combined/shRNA ]; then
		mkdir combined/shRNA
	fi
	if [[ -e combined/shRNA/tmp_shRNA_clusters_${analysisname}.bed ]]; then
		rm -f combined/shRNA/tmp_shRNA_clusters_${analysisname}.bed
	fi
	for tissue in ${uniq_shrna_tissue_list[@]}
	do
		printf "\nMaking shRNA cluster file for ${tissue}\n"
		awk -v OFS="\t" '($0 !~ /^#/) && ($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/) {print $1,$4,$5,$3"_"NR,"."}' shRNA/mapped/${line}_${tissue}_shRNA/ShortStack_All.gff3 | bedtools sort -g ${ref_dir}/chrom.sizes > combined/shRNA/${line}_${tissue}_shRNA_clusters.bed
		awk -v OFS="\t" -v t=${tissue} '{print $1,$2,$3,t}' combined/shRNA/${line}_${tissue}_shRNA_clusters.bed | sort -k1,1 -k2,2n -u >> combined/shRNA/tmp_shRNA_clusters_${analysisname}.bed
		printf "\nGetting closest gene for ${tissue}\n"	
		bedtools closest -a combined/shRNA/${line}_${tissue}_shRNA_clusters.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/shRNA/temp_clusters_${analysisname}.bed
		printf "\nGrouping based on distance for ${tissue}\n"
		awk -v OFS="\t" '{if ($5<-2000) {d="Intergenic"} else if ($5<0) {d="Terminator"} else if ($5==0) {d="Gene_body"} else if ($5>2000) {d="Intergenic"} else {d="Promoter"} print $0,d}' combined/shRNA/temp_clusters_${analysisname}.bed > combined/shRNA/temp2_clusters_${analysisname}.bed
		printf "\nIntersecting TE for ${tissue}\n"
		bedtools intersect -a combined/shRNA/temp2_clusters_${analysisname}.bed -b combined/TSS/${ref}_all_tes.bed -loj | awk -v OFS="\t" -v t=${tissue} -v l=${line} '{if ($13==".") print l,t,l"_"t"_"$4,$9,"No",$9,$9; else if ($9 == "Intergenic") print l,t,l"_"t"_"$4,$9,$13,$13,$13; else print l,t,l"_"t"_"$4,$9,$13,$13,$13"_in_"$9}' > combined/shRNA/Clusters_in_genes_and_tes_${tissue}_${analysisname}.bed
		rm -f combined/shRNA/temp*_clusters_${analysisname}.bed
	done
	printf "Line\tTissue\tCluster_ID\tGene\tTE\tLabel\tLabelcombined\n" > combined/shRNA/Table_shRNA_clusters_tissues_${analysisname}.txt
	cat combined/shRNA/Clusters_in_genes_and_tes_*_${analysisname}.bed >> combined/shRNA/Table_shRNA_clusters_tissues_${analysisname}.txt
	printf "\nPreparing merged shRNA cluster file for ${analysisname}\n"
	sort -k1,1 -k2,2n combined/shRNA/tmp_shRNA_clusters_${analysisname}.bed > combined/shRNA/tmp2_shRNA_clusters_${analysisname}.bed
	bedtools merge -i combined/shRNA/tmp2_shRNA_clusters_${analysisname}.bed -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{print $1,$2,$3,"Cluster_"NR,$4}'> combined/shRNA/tmp3_shRNA_clusters_${analysisname}.bed
	printf "\nGetting closest gene for clusters in ${analysisname}\n"
	if [[ ${ref} == "B73_v4" ]]; then
		bedtools closest -a combined/shRNA/tmp3_shRNA_clusters_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/shRNA/tmp4_shRNA_clusters_${analysisname}.bed
	else
		bedtools closest -a combined/shRNA/tmp3_shRNA_clusters_${analysisname}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > combined/shRNA/tmp4_shRNA_clusters_${analysisname}.bed
	fi
	printf "\nGrouping based on distance\n"
	awk -v OFS="\t" '{if ($5<-2000) {d="Intergenic"} else if ($5<0) {d="Terminator"} else if ($5==0) {d="Gene_body"} else if ($5>2000) {d="Intergenic"} else {d="Promoter"} print $0,d}' combined/shRNA/tmp4_shRNA_clusters_${analysisname}.bed > combined/shRNA/tmp5_shRNA_clusters_${analysisname}.bed
	printf "\nIntersecting with TEs\n"
	bedtools intersect -a combined/shRNA/tmp5_shRNA_clusters_${analysisname}.bed -b combined/TSS/${ref}_all_tes.bed -loj | awk -v OFS="\t" -v l=${line} 'BEGIN {printf "Line\tCluster_ID\tGene\tTE\tLabel\tLabelcombined\tGID\tTissues\tPeak_coordinates\n"} {if ($13==".") print l,$4,$9,"No",$9,$9,$8,$7,$1":"$2"-"$3; else if ($9 == "Intergenic") print l,$4,$9,$13,$13,$13,$8,$7,$1":"$2"-"$3; else print l,$4,$9,$13,$13,$13"_in_"$9,$8,$7,$1":"$2"-"$3}' > combined/shRNA/all_shRNA_clusters_in_genes_and_tes_${analysisname}.bed
	rm -f combined/shRNA/tmp*_shRNA_clusters_${analysisname}.bed
	#### To create a matrix of peak presence in each sample
	printf "\nCreating matrix file for ${analysisname}\n"
	for tissue in ${uniq_shrna_tissue_list[@]}
	do
		printf "${tissue}\n" > combined/shRNA/temp_col_clusters_${analysisname}_${tissue}.txt
		awk -v OFS="\t" -v t=${tissue} 'NR>1 {if ($8 ~ t) print "1"; else print "0"}' combined/shRNA/all_shRNA_clusters_in_genes_and_tes_${analysisname}.bed >> combined/shRNA/temp_col_clusters_${analysisname}_${tissue}.txt
	done
	awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7}' combined/shRNA/all_shRNA_clusters_in_genes_and_tes_${analysisname}.bed > combined/shRNA/temp_col_clusters_${analysisname}_AAA.txt
	paste combined/shRNA/temp_col_clusters_${analysisname}_*.txt | uniq > combined/shRNA/matrix_upset_shRNA_clusters_${analysisname}.txt
	rm -f combined/shRNA/temp_col_clusters_${analysisname}_*.txt
	#### To make an Upset plot highlighting peaks in gene bodies
	printf "\nCreating Distirbution and Upset plot for shRNA clusters in ${analysisname} with R version:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_shRNA_distribution_upset.r ${analysisname} ${TEtypestring} combined/shRNA/Table_shRNA_clusters_tissues_${analysisname}.txt combined/shRNA/matrix_upset_shRNA_clusters_${analysisname}.txt combined/shRNA/all_shRNA_clusters_in_genes_and_tes_${analysisname}.bed combined/DEG/genes_rpkm_${analysisname}.txt
fi

#########################################################################################
####################################### PART16 ##########################################
########################## Making heatmaps on all TEs ###################################
#########################################################################################

rnaseqsamples=${#rnaseq_bw_list_plus[@]}
rampagesamples=${#rampage_bw_list_plus[@]}
shrnasamples=${#shrna_bw_list_plus[@]}
totsamples=$((rnaseqsamples+rampagesamples+shrnasamples))
if [[ ${tefilebw} != "" ]] && [ ${totsamples} -gt 0 ] && [[ "${repeats}" == "YES" ]]; then
	#### Computing the stranded matrix
	while read TEtype
	do
		awk -v t=${TEtype} '$4==t && $6=="+"' combined/TSS/${ref}_all_tes.bed > combined/TSS/${ref}_${TEtype}_${analysisname}_plus.bed
		awk -v t=${TEtype} '$4==t && $6=="-"' combined/TSS/${ref}_all_tes.bed > combined/TSS/${ref}_${TEtype}_${analysisname}_minus.bed
		for strand in plus minus
		do
			case "${strand}" in
				plus) 	bw_list="${rnaseq_bw_list_plus[@]} ${rampage_bw_list_plus[@]} ${shrna_bw_list_plus[@]}";;
				minus) 	bw_list="${rnaseq_bw_list_minus[@]} ${rampage_bw_list_minus[@]} ${shrna_bw_list_minus[@]}";;
			esac
			nb=$(wc -l combined/TSS/${ref}_${TEtype}_${analysisname}_${strand}.bed | awk '{print $1}')
			if [[ ${nb} -gt 0 ]]; then
				printf "\nComputing scale-regions ${strand} strand matrix for ${TEtype} from ${analysisname}\n"
				computeMatrix scale-regions -q --missingDataAsZero --skipZeros -R combined/TSS/${ref}_${TEtype}_${analysisname}_${strand}.bed -S ${bw_list} -bs 100 -b 2000 -a 2000 -m 5000 -p ${threads} -o combined/matrix/temp_TE_regions_${TEtype}_${analysisname}_${strand}.gz
				printf "\nComputing reference-point on TSS ${strand} strand matrix for ${TEtype} from ${analysisname}\n"
				computeMatrix reference-point --referencePoint "TSS" -q --missingDataAsZero --skipZeros -R combined/TSS/${ref}_${TEtype}_${analysisname}_${strand}.bed -S ${bw_list} -bs 100 -b 2000 -a 8000 -p ${threads} -o combined/matrix/temp_TE_tss_${TEtype}_${analysisname}_${strand}.gz
			fi
		done
		#### Merging stranded matrix, extracting scales and plotting heatmaps
		all_samples=()
		all_labels=()
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
		for matrix in TE_regions TE_tss
		do
			mat="empty"
			if [[ -s combined/matrix/${matrix}_${TEtype}_${analysisname}_plus.gz ]] && [[ -s combined/matrix/${matrix}_${TEtype}_${analysisname}_minus.gz ]]; then
				printf "\nMerging stranded matrices aligned by ${matrix} ${TEtype} of ${analysisname}\n"
				computeMatrixOperations rbind -m combined/matrix/${matrix}_${TEtype}_${analysisname}_plus.gz combined/matrix/${matrix}_${TEtype}_${analysisname}_minus.gz -o combined/matrix/${matrix}_${TEtype}_${analysisname}.gz
				mat="combined/matrix/${matrix}_${TEtype}_${analysisname}.gz"
			elif [[ -s combined/matrix/${matrix}_${TEtype}_${analysisname}_plus.gz ]]; then
				mat="combined/matrix/${matrix}_${TEtype}_${analysisname}_plus.gz"
			elif [[ -s combined/matrix/${matrix}_${TEtype}_${analysisname}_minus.gz ]]; then
				mat="combined/matrix/${matrix}_${TEtype}_${analysisname}_minus.gz"
			fi
			if [[ ${mat} != "empty" ]]; then	
				printf "\nGetting scales for ${matrix} ${TEtype} matrix of ${analysisname}\n"
				computeMatrixOperations dataRange -m ${mat} > combined/matrix/temp_values_${matrix}_${TEtype}_${analysisname}.txt
				plotProfile -m ${mat} -out combined/plots/temp_${matrix}_${TEtype}_${analysisname}_profile.pdf --samplesLabel ${all_labels[@]} --averageType mean --outFileNameData combined/matrix/temp_values_profile_${matrix}_${TEtype}_${analysisname}.txt
				rm -f combined/plots/temp_${matrix}_${TEtype}_${analysisname}_profile.pdf
				mins=()
				maxs=()
				ymins=()
				ymaxs=()
				for mark in ${all_samples[@]}
				do
					mini=$(grep "${mark}" combined/matrix/temp_values_${matrix}_${TEtype}_${analysisname}.txt | awk 'BEGIN {m=999999} {a=$5; if (a<m) m=a;} END {print m}')
					maxi=$(grep "${mark}" combined/matrix/temp_values_${matrix}_${TEtype}_${analysisname}.txt | awk 'BEGIN {m=-999999} {a=$6; if (a>m) m=a;} END {print m}')
					num=$(grep "${mark}" combined/matrix/temp_values_${matrix}_${TEtype}_${analysisname}.txt | wc -l)
					test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
					if [[ ${test} == "yes" ]]; then
						mini=("0")
						maxi=("0.005")
					fi
					for i in $(seq 1 ${num})
					do
						mins+=("${mini}")
						maxs+=("${maxi}")
					done		
					ymini=$(grep "${mark}" combined/matrix/temp_values_profile_${matrix}_${TEtype}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
					ymaxi=$(grep "${mark}" combined/matrix/temp_values_profile_${matrix}_${TEtype}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
					num=$(grep "${mark}" combined/matrix/temp_values_profile_${matrix}_${TEtype}_${analysisname}.txt | wc -l)
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
					mini=$(grep ${sample} combined/matrix/temp_values_${matrix}_${TEtype}_${analysisname}.txt | awk '{print $5}')
					maxi=$(grep ${sample} combined/matrix/temp_values_${matrix}_${TEtype}_${analysisname}.txt | awk '{print $6}')
					test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
					if [[ ${test} == "yes" ]]; then
						mins2+=("0")
						maxs2+=("0.005")
					else
						mins2+=("${mini}")
						maxs2+=("${maxi}")
					fi
				done
				ymins2=()
				ymaxs2=()
				for sample in ${all_labels[@]}
				do
					ymini=$(grep ${sample} combined/matrix/temp_values_profile_${matrix}_${TEtype}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
					ymaxi=$(grep ${sample} combined/matrix/temp_values_profile_${matrix}_${TEtype}_${analysisname}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
					test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
					if [[ ${test} == "yes" ]]; then
						ymins2+=("0")
						ymaxs2+=("0.01")
					else
						ymins2+=("${ymini}")
						ymaxs2+=("${ymaxi}")
					fi
				done
				printf "\nPlotting heatmap for ${matrix} ${TEtype} matrix of ${analysisname} scaling by mark\n"
				plotHeatmap -m ${mat} -out combined/plots/${analysisname}_heatmap_${matrix}_${TEtype}.pdf --sortRegions descend --sortUsing mean --samplesLabel ${all_labels[@]} --colorMap 'seismic' --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear'
				printf "\nPlotting heatmap for ${matrix} ${TEtype} matrix of ${analysisname} scaling by sample\n"
				plotHeatmap -m ${mat} -out combined/plots/${analysisname}_heatmap_${matrix}_${TEtype}_v2.pdf --sortRegions descend --sortUsing mean --samplesLabel ${all_labels[@]} --colorMap 'seismic' --zMin ${mins2[@]} --zMax ${maxs2[@]} --yMin ${ymins2[@]} --yMax ${ymaxs2[@]} --interpolationMethod 'bilinear'
			fi
		done
	done < combined/TSS/${ref}_TE_types.txt
	rm -f combined/matrix/temp*${analysisname}*
fi

####################

printf "\nCombined analysis script finished successfully for ${analysisname}\n"
touch combined/chkpts/analysis_${analysisname}
