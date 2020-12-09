#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=0.5G
#$ -l tmp_free=1G
#$ -o maizecodeanalysis.log
#$ -j y
#$ -N maizecodeanalysis

usage="
##### Script for Maize code data analysis
#####
##### sh MaiCode_analysis.sh -f samplefile [-r regionfile] [-s]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Sample, PE or SE, Reference genome directory
##### 	-r: textfile containing the name of region files that are to be plotted over (bed files)
##### 		It is safest to use a full paths.
#####		If no region file is given, the analysis will behave as if -s was set.
#####	-s: If set, the script does not progress into the line data analysis, only single sample analysis will be performed
##### 	-h: help, returns usage
##### 
##### It sends each type of sample to its specific analysis file (MaizeCode_ChIP_analysis.sh or MaizeCode_RNA_analysis.sh)
##### Then starts an analysis per reference line (MaizeCode_line_analysis.sh)
#####
##### Requirements: samtools, bedtools, deeptools, macs2, idr, R (+R packages: ggplot2,UpSetR,limma,edgeR,dplyr,tidyr,stringr,gplots,RColorBrewer,cowplot)
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS
export mc_dir="${PWD}/scripts/"
printf "\nRunning MaizeCode_analysis.sh script in working directory ${PWD}\n"

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
		*)	printf "\nArgument unknown, retunring usage:\n$usage\n"
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

if [ -s combined/temp_reports_${samplename}_ChIP.txt ]; then
	rm -f combined/temp_reports_${samplename}_ChIP.txt
fi

if [ -s combined/temp_reports_${samplename}_RNA.txt ]; then
	rm -f combined/temp_reports_${samplename}_RNA.txt
fi

#### Check if there are new samples to analyze individually
new_chipsample=()
new_rnasample=()
datatype_list=()
ref_list=()
while read line tissue sample paired ref_dir
do
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
		*) datatype="unknown";;
	esac
	name=${line}_${tissue}_${sample}
	printf "$line\t$tissue\t$sample\t$paired\t${ref_dir}\n" >> combined/temp_reports_${samplename}_${datatype}.txt
	if [ -e $datatype/chkpts/analysis_${name} ]; then
		printf "\nSingle sample analysis for $name already done!\n"	
	elif [[ "$datatype" == "ChIP" ]]; then
		if [ ! -d ./ChIP/peaks ]; then
			mkdir ./ChIP/peaks
		fi
		if [ ! -d ./ChIP/plots ]; then
			mkdir ./ChIP/plots
		fi
		datatype_list+=("${datatype}")
		new_chipsample+=("${name}")
		printf "$line\t$tissue\t$sample\t$paired\t${ref_dir}\n" >> $datatype/temp_${samplename}_${datatype}.txt
	elif [[ "$datatype" == "RNA" ]]; then
		if [ ! -d ./RNA/plots ]; then
			mkdir ./RNA/plots
		fi
		if [ ! -d ./RNA/TSS ]; then
			mkdir ./RNA/TSS
		fi
		datatype_list+=("${datatype}")
		new_rnasample+=("${name}")
		printf "$line\t$tissue\t$sample\t$paired\t${ref_dir}\n" >> $datatype/temp_${samplename}_${datatype}.txt
	else
		printf "\nType of data unknown for $name\nSample not processed\n"
	fi
	ref=${ref_dir##*/}
	if [[ ! "${ref_list[@]}" =~ "${ref}" ]]; then
		ref_list+=("$ref")
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
	wait ${pids[*]}

	for datatype in ${uniq_datatype_list[@]}
	do
		if [[ "$datatype" == "ChIP" ]]; then
			for chipsample in ${new_chipsample[@]}
			do
				if [ ! -e ${datatype}/chkpts/analysis_${chipsample} ]; then
					printf "\nProblem during the processing of ChIP sample ${chipsample}!\nCheck log: ChIP/logs/${samplename}.log and ChIP/logs/analysis_${chipsample}_*_.log\n"
				else
					printf "\nChIP analysis for $chipsample processed succesfully\n"
				fi
			done
		elif [[ "$datatype" == "RNA" ]]; then
			for rnasample in ${new_rnasample[@]}
			do
				if [ ! -e ${datatype}/chkpts/analysis_${rnasample} ]; then
					printf "\nProblem during the processing of RNA sample ${rnasample}!\nCheck log: RNA/logs/${samplename}.log and RNA/logs/analysis_${rnasample}.log\n"
				else 
					printf "\nRNA analysis for $rnasample processed succesfully\n"
				fi
			done
		fi
	done
fi

#### To get the peaks stats for all ChIPseq samples in the samplefile

if [ -s combined/temp_reports_${samplename}_ChIP.txt ]; then
	printf "\nSummarizing peak stats for ${samplename}\n"
	if [ -s combined/reports/temp_peaks_${samplename}.txt ]; then
		rm -f combined/reports/temp_peaks_${samplename}.txt
	fi
	while read line tissue mark paired ref_dir
	do
		awk -v a=$line -v b=$tissue -v c=$mark '$1==a && $2==b && $3==c' ChIP/reports/summary_ChIP_peaks.txt >> combined/reports/temp_peaks_${samplename}.txt
	done < combined/temp_reports_${samplename}_ChIP.txt
	printf "Line\tTissue\tMark\tPeaks_in_Rep1\tPeaks_in_Rep2\tCommon_peaks\tCommon_peaks_IDR_0.05\tPeaks_in_merged\tPeaks_in_pseudo_reps\tSelected_peaks\n" > combined/reports/summary_ChIP_peaks_${samplename}.txt
	sort combined/reports/temp_peaks_${samplename}.txt -u >> combined/reports/summary_ChIP_peaks_${samplename}.txt
	rm -f combined/reports/temp_peaks_${samplename}.txt
	printf "\nPlotting peak stats for all samples in the samplefile with R:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_peak_stats.r combined/reports/summary_ChIP_peaks_${samplename}.txt ${samplename}
fi

#### To get the RNA stats for all RNA samples in the samplefile

if [ -s combined/temp_reports_${samplename}_RNA.txt ]; then
	#### To get gene expression stats for RNAseq samples
	grep "RNAseq" combined/temp_reports_${samplename}_RNA.txt > combined/reports/temp_${samplename}.txt
	exist=$(cat combined/reports/temp_${samplename}.txt | wc -l)
	if [ $exist -gt 0 ]; then
		printf "\nSummarizing gene expression stats for ${samplename}\n"
		if [ -s combined/reports/temp_gene_expression_${samplename}.txt ]; then
			rm -f combined/reports/temp_gene_expression_${samplename}.txt
		fi
		while read line tissue sample paired ref_dir
		do
			awk -v a=$line -v b=$tissue -v c=$sample '$1==a && $2==b && $3==c' RNA/reports/summary_gene_expression.txt >> combined/reports/temp_gene_expression_${samplename}.txt
		done < combined/reports/temp_${samplename}.txt
		printf "Line\tTissue\tType\tTotal_annotated_genes\tNot_expressed_in_Rep1\tLow_expression_in_Rep1(<1cpm)\tHigh_expression_in_Rep1(>1cpm)\tNot_expressed_in_Rep2\tLow_expression_in_Rep2(<1cpm)\tHigh_expression_in_Rep2(>1cpm)\tNo_mean\tLow_mean(<1cpm)\tHigh_mean(>1cpm)\n" > combined/reports/summary_gene_expression_${samplename}.txt
		sort combined/reports/temp_gene_expression_${samplename}.txt -u >> combined/reports/summary_gene_expression_${samplename}.txt
		rm -f combined/reports/temp_gene_expression_${samplename}.txt
		printf "\nPlotting gene expression stats for all RNAseq samples in the samplefile with R:\n"
		R --version
		Rscript --vanilla ${mc_dir}/MaizeCode_R_gene_ex_stats.r combined/reports/summary_gene_expression_${samplename}.txt ${samplename}
	fi
	#### To get tss stats for RAMPAGE samples
	grep "RAMPAGE" combined/temp_reports_${samplename}_RNA.txt > combined/reports/temp_${samplename}.txt
	exist=$(cat combined/reports/temp_${samplename}.txt | wc -l)
	if [ $exist -gt 0 ]; then
		printf "\nSummarizing tss stats for ${samplename}\n"
		if [ -s combined/reports/temp_RAMPAGE_tss_${samplename}.txt ]; then
			rm -f combined/reports/temp_RAMPAGE_tss_${samplename}.txt
		fi
		while read line tissue sample paired ref_dir
		do
			awk -v a=$line -v b=$tissue -v c=$sample '$1==a && $2==b && $3==c' RNA/reports/summary_RAMPAGE_tss.txt >> combined/reports/temp_RAMPAGE_tss_${samplename}.txt
		done < combined/reports/temp_${samplename}.txt
		printf "Line\tTissue\tType\tTotal_annotated_genes\tTSS_in_rep1\tTSS_in_Rep2\tCommon_TSS\tCommon_TSS_IDR<=0.05\n" > combined/reports/summary_RAMPAGE_tss_${samplename}.txt
		sort combined/reports/temp_RAMPAGE_tss_${samplename}.txt -u >> combined/reports/summary_RAMPAGE_tss_${samplename}.txt
		rm -f combined/reports/temp_RAMPAGE_tss_${samplename}.txt
	fi
	if [ -e combined/reports/temp_${samplename}.txt ]; then
		rm -f combined/reports/temp_${samplename}.txt
	fi
fi

rm -f combined/temp_reports_*

if [[ $keepgoing == "STOP" ]]; then
	printf "\nScript finished successfully without combined analysis\n"
	touch chkpts/${analysisname}
	exit 0
fi	

#############################################################################################
########################################### PART2 ###########################################
#################### Do the analysis of ChIP and RNA data for each reference ################
#############################################################################################

if [ -e combined/chkpts/${analysisname} ]; then
	rm -f combined/chkpts/${analysisname}
fi

check_list=()
region_list=()
pids=()
##### Processing line analysis per reference genome
for ref in ${ref_list[@]}
do
	grep "$ref" $samplefile > combined/${samplename}_analysis_samplefile.temp_${ref}.txt
	grep "$ref" $regionfile > combined/${regionname}.temp_${ref}.txt
	regioni=$(cat combined/${regionname}.temp_${ref}.txt)
	tmp3=${regioni##*/}
	regioniname=${tmp3%%.*}
	check_list+=("combined/chkpts/analysis_${samplename}_on_${regioniname}")
	region_list+=("${regioniname}")
	printf "\nLaunching line analysis script for samplefile $samplename on regionfile $regioniname\n"
	qsub -sync y -N ${ref}_analysis -o combined/logs/analysis_${samplename}_on_${regioniname}.log ${mc_dir}/MaizeCode_line_analysis.sh -f combined/${samplename}_analysis_samplefile.temp_${ref}.txt -r ${regioni} &
	pids+=("$!")
done

printf "\nWaiting for line analysis to be processed\n"
wait ${pids[*]}

i=0
for check in ${check_list[@]}
do
	if [ ! -e ${check} ]; then
		printf "\nProblem during the line analysis of ${ref_list[i]} on ${region_list[i]}!\nCheck log: combined/logs/analysis_${samplename}_on_${region_list[i]}.log\n"
	else 
		printf "\nLine analysis of ${ref_list[i]} on ${region_list[i]} processed succesfully\n"
		rm -f combined/${samplename}_analysis_samplefile.temp_${ref_list[i]}.txt
		rm -f combined/${regionname}.temp_${ref_list[i]}.txt
	fi
	i=$((i+1))
done

rm -f combined/*temp*
touch combined/chkpts/${analysisname}


#############################################################################################
########################################### PART3 ###########################################
#################### Do the analysis of ChIP and RNA data accross references ################
#############################################################################################


#############################################################################################
########################################### MISC ############################################
#############################################################################################

###### To make the samplefile for B73_endosperm

# printf "B73\tendosperm\tH3K4me1\tPE\t/grid/martienssen/home/jcahn/data/Genomes/Zea_mays/B73_v4\nB73\tendosperm\tH3K4me3\tPE\t/grid/martienssen/home/jcahn/data/Genomes/Zea_mays/B73_v4\nB73\tendosperm\tH3K27ac\tPE\t/grid/martienssen/home/jcahn/data/Genomes/Zea_mays/B73_v4\nB73\tendosperm\tRNAseq\tPE\t/grid/martienssen/home/jcahn/data/Genomes/Zea_mays/B73_v4\nB73\tendosperm\tRAMPAGE\tPE\t/grid/martienssen/home/jcahn/data/Genomes/Zea_mays/B73_v4" > B73_endosperm_analysis_samplefile.txt

#############################################################################################

