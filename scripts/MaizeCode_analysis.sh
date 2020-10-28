#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=2G
#$ -l tmp_free=10G
#$ -o maizecodeanalysis.log
#$ -j y
#$ -N maizecodeanalysis

usage="
##### Script for Maize code data analysis
#####
##### sh MaiCode_analysis.sh -f samplefile [-r regionfile] [-s]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Sample, Rep (Rep1, Rep2 or merged), PE or SE
##### 	-r: bedfile containing the regions that want to be plotted over
##### 		(safest to use a full path to the region file)
#####		If no region file is given, the analysis will behave as if -s was set
#####	-s: If set, the script does not progress into the combined data analysis
##### 	-h: help, returns usage
##### 
##### It sends each type of sample to its specific analysis file (MaizeCode_ChIP_analysis.sh or MaizeCode_RNA_analysis.sh)
##### Then starts a combined analysis (MaizeCode_combined_analysis.sh)
##### If replicates are to be merged, 'merged' should be the value in column #4 (replicate ID) of the samplefile
##### For cleaner naming purposes, use 'analysis_samplefile.txt' as suffix
#####
##### Requirements: samtools, bedtools, deeptools, macs2, idr, R (+R packages: ggplot2,readr,UpSetR)
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS
# # export mc_dir=$(dirname "$0")
export mc_dir="${HOME}/data/Scripts/MaizeCode"
printf "\nRunning MaizeCode scripts from ${mc_dir} in working directory ${PWD}\n"

if [ $# -eq 0 ]; then
	printf "$usage\n"
	exit 1
fi

while getopts ":f:r:sh" opt; do
	case $opt in
		f) 	export samplefile=${OPTARG};;
		r)	export regionfile=${OPTARG};;
		s)	printf "\nOption not to perform combined analysis selected\n"
			export keepgoing="STOP";;
		h) 	printf "$usage\n"
			exit 0;;
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
	printf "Regionfile is missing.\nAnalysis will be stopped without deeper analysis\nIf that was not intended, check usage.\n"
	keepgoing="STOP"
fi

#############################################################################################
########################################### PART1 ###########################################
######################## Send each sample to specific analysis type #########################
#############################################################################################

tmp1=${samplefile##*/}
samplename=${tmp1%%_analysis*}

if [[ "$keepgoing" == "STOP" ]]; then
	analysisname="${samplename}_no_region"
else
	tmp2=${regionfile##*/}
	regionname=${tmp2%.*}
	analysisname="${samplename}_on_${regionname}"
fi

if [ -s ChIP/temp_${samplename}_ChIP.txt ]; then
	rm -f ChIP/temp_${samplename}_ChIP.txt
fi

if [ -s RNA/temp_${samplename}_RNA.txt ]; then
	rm -f RNA/temp_${samplename}_RNA.txt
fi

#### Check if there are new samples to analyze individually
new_chipsample=()
new_rnasample=()
datatype_list=()
while read line tissue sample rep paired
do
	name=${line}_${tissue}_${sample}_${rep}
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
		*) datatype="unknown";;
	esac
	if [ -e $datatype/chkpts/analysis_${name} ]; then
		printf "\nSingle sample analysis for $name already done!\n"	
	elif [[ "$datatype" == "ChIP" ]]; then
		datatype_list+=("${datatype}")
		new_chipsample+=("${name}")
		if [ ! -d ./ChIP/peaks ]; then
			mkdir ./ChIP/peaks
		fi
		if [ ! -d ./ChIP/plots ]; then
			mkdir ./ChIP/plots
		fi
		printf "$line\t$tissue\t$sample\t$rep\t$paired\n" >> $datatype/temp_${samplename}_ChIP.txt
	elif [[ "$datatype" == "RNA" ]]; then
		datatype_list+=("${datatype}")
		new_rnasample+=("${name}")
		printf "$line\t$tissue\t$sample\t$rep\t$paired\n" >> $datatype/temp_${samplename}_RNA.txt
	else
		printf "\nType of data unknown for $name\nSample not processed\n"
	fi
done < $samplefile

#### If there are new samples, run the ChIP and/or RNA pipeline on the new samples of the same type

test_new=1
if [ ${#new_chipsample[@]} -eq 0 ] && [ ${#new_rnasample[@]} -eq 0 ]; then
	printf "\nAll samples in samplefile already analyzed individually\n"
	test_new=0
fi

if [[ "${test_new}" == 1 ]]; then
	uniq_datatype_list=($(printf "%s\n" "${datatype_list[@]}" | sort -u))
	pids=()
	for datatype in ${uniq_datatype_list[@]}
	do
		printf "\nRunning $datatype analysis script\n"
		cd $datatype
		qsub -sync y -N ${datatype} -o logs/${samplename}.log ${mc_dir}/MaizeCode_${datatype}_analysis.sh -f temp_${samplename}_${datatype}.txt &
		pids+=("$!")
		cd ..
	done
	#### Wait for the other scripts to finish
	printf "\nWaiting for the datatype analysis scripts to finish\n"
	wait $pids

	for datatype in ${uniq_datatype_list[@]}
	do
		if [[ "$datatype" == "ChIP" ]]; then
			for chipsample in ${new_chipsample[@]}
			do
				if [ ! -e ${datatype}/chkpts/analysis_${chipsample} ]; then
					printf "\nProblem during the processing of ChIP sample ${chipsample}!\nCheck log: ChIP/logs/${samplename}.log\n"
				else 
					printf "\nChIP analysis for $chipsample processed succesfully\n"
				fi
			done
		elif [[ "$datatype" == "RNA" ]]; then
			for rnasample in ${new_rnasample[@]}
			do
				if [ ! -e ${datatype}/chkpts/analysis_${rnasample} ]; then
					printf "\nProblem during the processing of RNA sample ${rnasample}!\nCheck log: RNA/logs/${samplename}.log\n"
				else 
					printf "\nChIP analysis for $rnasample processed succesfully\n"
				fi
			done
		fi
	done
fi

if [[ $keepgoing == "STOP" ]]; then
	printf "\nScript finished successfully without combined analysis\n"
	touch chkpts/${analysisname}
	exit 0
fi	

#############################################################################################
########################################### PART2 ###########################################
###################### Do the combined analysis of ChIP and RNA data ########################
#############################################################################################

if [ ! -d ./combined ]; then
	mkdir ./combined
	mkdir ./combined/peaks
	mkdir ./combined/DEG
	mkdir ./combined/matrix
	mkdir ./combined/plots
	mkdir ./combined/chkpts
	mkdir ./combined/logs
fi

if [ -e combined/chkpts/combined_${analysisname} ]; then
	rm -f combined/chkpts/${analysisname}
fi

pids=()
##### Processing combined analysis
printf "\nLaunching combined analysis script\n"
qsub -sync y -N combined_analysis -o combined/logs/combined_analysis.log ${mc_dir}/MaizeCode_combined_analysis.sh -f $samplefile -r $regionfile &
pids+=("$!")

wait $pids

if [ ! -e combined/chkpts/${analysisname} ]; then
	printf "\nProblem during the combined analysis of ${samplename} on ${regionname}!\nCheck log: logs/${analysisname}.log\n"
else 
	printf "\nCombined analysis for $samplename on $regionname processed succesfully\n"
	touch combined/chkpts/${analysisname}
fi

#############################################################################################
########################################### MISC ############################################
#############################################################################################

###### To make a test samplefile for (B73 endosperm H3K4me1)

# printf "B73\tendosperm\tH3K4me1\tRep1\tPE\nB73\tendosperm\tH3K4me1\tRep2\tPE\n" > test_analysis_samplefile.txt

###### To make the samplefile for B73_endosperm

# printf "B73\tendosperm\tH3K4me1\tRep1\tPE\nB73\tendosperm\tH3K4me1\tRep2\tPE\nB73\tendosperm\tH3K4me1\tmerged\tPE\nB73\tendosperm\tH3K4me3\tRep1\tPE\nB73\tendosperm\tH3K4me3\tRep2\tPE\nB73\tendosperm\tH3K4me3\tmerged\tPE\nB73\tendosperm\tH3K27ac\tRep1\tPE\nB73\tendosperm\tH3K27ac\tRep2\tPE\nB73\tendosperm\tH3K27ac\tmerged\tPE\n" > B73_endosperm_analysis_samplefile.txt

###### To create a regionfile containing several groups of regions (e.g. gene_list1, gene_list2 and gene_list3)

# printf "#Name_gene_list1" > regionfile.bed
# cat gene_list1.bed >> regionfile.bed
# printf "#Name_gene_list2" >> regionfile.bed
# cat gene_list2.bed >> regionfile.bed
# printf "#Name_gene_list3" >> regionfile.bed
# cat gene_list3.bed >> regionfile.bed

#############################################################################################

