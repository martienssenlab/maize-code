#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 2
#$ -l m_mem_free=2G
#$ -l tmp_free=10G
#$ -o shRNAanalysis.log
#$ -j y
#$ -N shRNAanalysis

usage="
##### Script for Maize code shRNA data analysis, used by script MaizeCode_analysis.sh for shRNA data
#####
##### sh MaiCode_shRNA_analysis.sh -f samplefile [-h]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Mark, PE or SE, Reference genome directory
##### 	-h: help, returns usage
##### 
##### It merges the two replicate files (filtered for sizes from 15 to 32 nt), maps it again with ShortStack and creates stranded bigwig files for each sample
#####
##### Requirements: samtools, deeptools, Shortstack
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS

if [ $# -eq 0 ]; then
	printf "${usage}\n"
	exit 1
fi

while getopts ":f:h" opt; do
	case ${opt} in
		h) 	printf "${usage}\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		*)	printf "${usage}\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! ${samplefile} ]; then
	printf "Samplefile missing!\n"
	printf "${usage}\n"
	exit 1
fi

tmp1=${samplefile##*temp_}
export samplename=${tmp1%%_shRNA}

pids=()
while read line tissue rnatype paired ref_dir
do
	export line
	export tissue
	export rnatype
	export name=${line}_${tissue}_${rnatype}
	printf "\nStarting single RNA sample analysis for ${name}\n"	
	export ref_dir=${ref_dir}
	export ref=${ref_dir##*/}
	if [ -s ${ref_dir}/*.fa ]; then
		fa_file=$(ls ${ref_dir}/*.fa)
		fasta=${fa_file}
	elif [ -s ${ref_dir}/*.fasta ]; then
		fa_file=$(ls ${ref_dir}/*.fasta)
		fasta=${fa_file}
	elif [ -s ${ref_dir}/*.fa.gz ]; then
		fa_file=$(ls ${ref_dir}/*.fa.gz)
		pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${rnatype}_${ref}.fa
		fasta=${ref_dir}/temp_${rnatype}_${ref}.fa
	elif [ -s ${ref_dir}/*.fasta.gz ]; then
		fa_file=$(ls ${ref_dir}/*.fasta.gz)
		pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${rnatype}_${ref}.fa
		fasta=${ref_dir}/temp_${rnatype}_${ref}.fa	
	fi
	export fasta
	qsub -N ${name} -V -cwd -sync y -pe threads 20 -l m_mem_free=5G -l tmp_free=50G -j y -o logs/analysis_${name}.log <<-'EOF' &
		#!/bin/bash
		set -e -o pipefail
		
		printf "\n\n"
		date
		printf "\n"
		
		export threads=$NSLOTS
		
		#### To merge bam files of replicates
		if [ ! -s mapped/${name}_merged.bam ]; then
			printf "\nMerging replicates of ${name}\n"
			samtools merge -@ ${threads} mapped/temp_${name}_merged.bam mapped/${name}_Rep*/sized_${name}_Rep*.bam
			samtools sort -@ ${threads} -o mapped/${name}_merged.bam mapped/temp_${name}_merged.bam
			rm -f mapped/temp_${name}_merged.bam
			samtools index -@ ${threads} mapped/${name}_merged.bam
		fi
		#### To map and identify sRNA loci with ShortStack
		if [ ! -s mapped/${name}/ShortStack_All.gff3 ]; then
			rm -fr mapped/${name}
			ShortStack --bamfile mapped/${name}_merged.bam --genomefile ${fasta} --sort_mem 6G --mmap u --dicermin 20 --dicermax 24 --bowtie_m 1000 --mismatches 1 --foldsize 1000 --pad 250 --outdir mapped/${name}
		fi
		#### To make bw files of merged samples if not already existing		
		if [ ! -s tracks/${name}_merged_minus.bw ]; then
			printf "\nMaking tracks for ${name}\n"
			bamCoverage --filterRNAstrand reverse -bs 1 -p ${threads} --normalizeUsing CPM -b mapped/${name}_merged.bam -o tracks/${name}_merged_plus.bw
			bamCoverage --filterRNAstrand forward -bs 1 -p ${threads} --normalizeUsing CPM -b mapped/${name}_merged.bam -o tracks/${name}_merged_minus.bw
		else
			printf "\nBigwig files for ${name} already exists\n"
		fi
		printf "\nAnalysis finished for ${name}\n"
		touch chkpts/analysis_${name}
	EOF
	pids+=("$!")
done < ${samplefile}

printf "\nWaiting for samples to be processed individually\n"
wait ${pids[*]}

rm -f ${ref_dir}/temp_${rnatype}_${ref}*

printf "\nScript finished successfully!\n"
