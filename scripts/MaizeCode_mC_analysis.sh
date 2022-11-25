#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -l tmp_free=2G
#$ -o mCanalysis.log
#$ -j y
#$ -N mCanalysis

usage="
##### Script for Maize code DNA methylation data analysis, used by script MaizeCode_analysis.sh for mC data
#####
##### sh MaiCode_mC_analysis.sh -f samplefile [-h]
#####	-f: samplefile containing the samples to compare and the reference directory in 6 tab-delimited columns:
##### 		Data, Line, Tissue, Mark, PE or SE, Reference directory
##### 	-h: help, returns usage
##### 
##### It merges the two replicate files (weighted average), and creates bigwig files
#####
##### Requirements: bedtools, bedGraphToBigWig
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

while getopts ":f:a:h" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		a)	export mapparam=${OPTARG};;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $samplefile ]; then
	printf "Samplefile missing!\n"
	printf "$usage\n"
	exit 1
fi

if [ ! ${mapparam} ]; then
	printf "No mapping option selected, using default\n"
	export mapparam="default"
elif [[ "${mapparam}" == "default" ]] || [[ "${mapparam}" == "Colcen" ]]; then
	printf "${mapparam} chosen as the mapping option\n"
else
	printf "Unknown mapping option selected\n"
	printf "${usage}\n"
	exit 1
fi

if [ ! -s reports/summary_mC_coverage.txt ]; then
	printf "Line\tTissue\tMark\tPeaks_in_Rep1\tPeaks_in_Rep2\tCommon_peaks\tCommon_peaks_IDR_0.05\tPeaks_in_merged\tPeaks_in_pseudo_reps\tSelected_peaks\n" > reports/summary_mC_coverage.txt
fi

pidsa=()
while read data line tissue mark paired ref_dir
do
	#### To merge bam files of replicates
	export line
	export tissue
	export mark
	export ref_dir
	export name=${line}_${tissue}_${mark}
	printf "\nStarting single ChIP sample analysis for $name\n"
	qsub -N ${name} -V -cwd -sync y -pe threads 12 -l m_mem_free=2G -l tmp_free=4G -j y -o logs/analysis_${name}.log <<-'EOF1' &
		#!/bin/bash
		set -e -o pipefail		
		export threads=$NSLOTS
		
		printf "\nMerging replicates of $name\n"
		printf "${sample}\n"
		zcat methylcall/${name}_Rep*.deduplicated.CX_report.txt.gz | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2-1,$2,$3,$4,$5,$6,$7}' > methylcall/temp_${name}.bed
		bedtools merge -d -1 -o distinct,sum,sum,distinct,distinct -c 4,5,6,7,8 -i methylcall/temp_${name}.bed > methylcall/temp2_${name}.bed
		cat methylcall/temp2_${name}.bed | awk -v OFS="\t" -v s=${name} '($5+$6)>0 {a=$5+$6; if ($7=="CHH") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CHH.bedGraph"; else if ($7=="CHG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CHG.bedGraph"; else if ($7=="CG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CG.bedGraph"}'
		rm -f methylcall/temp*${name}*
		for context in CG CHG CHH
		do
			printf "\nMaking bigwig files of ${context} context for ${name}\n"
			LC_COLLATE=C sort -k1,1 -k2,2n methylcall/${name}_${context}.bedGraph > methylcall/sorted_${name}_${context}.bedGraph
			bedGraphToBigWig methylcall/sorted_${name}_${context}.bedGraph ${ref_dir}/chrom.sizes methylcall/${name}_${context}.bw
		done
		rm -f methylcall/*${name}*.bedGraph
		touch chkpts/analysis_${name}
	EOF1
	pidsa+=("$!")
done < $samplefile

printf "\nWaiting for each sample to be processed individually\n"
wait ${pidsa[*]}
printf "\nScript finished successfully!\n"
