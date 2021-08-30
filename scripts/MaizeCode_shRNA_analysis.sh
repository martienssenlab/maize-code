#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=6G
#$ -l tmp_free=8G
#$ -o temp.log
#$ -j y
#$ -N temp

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS
export mc_dir="${PWD}/scripts/"
printf "\nRunning MaizeCode.sh script in working directory ${PWD}\n"

if [ $# -eq 0 ]; then
	printf "$usage\n"
	exit 1
fi

while getopts "d:l:t:m:p:h" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		d) 	export ref_dir=${OPTARG};;
		l)	export line=${OPTARG};;
		t)	export tissue=${OPTARG};;
		m)	export rnatype=${OPTARG};;
		p)	export paired=${OPTARG};;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $ref_dir ] || [ ! $line ] || [ ! $tissue ] || [ ! $rnatype ] || [ ! $paired ]; then
	printf "Missing arguments!\n"
	printf "$usage\n"
	exit 1
fi

export genomes=${ref_dir%/*}
export ref=${ref_dir##*/}

name=${line}_${tissue}_${rnatype}
### merging bam files
if [ ! -s mapped/${name}_merged.bam ]; then
	printf "\nMerging bam files\n"
	samtools merge -@ $threads mapped/temp_${name}_merged.bam mapped/${name}_Rep1/filtered_${name}_Rep1.bam mapped/${name}_Rep2/filtered_${name}_Rep2.bam 
	samtools sort -@ $threads -o mapped/${name}_merged.bam mapped/temp_${name}_merged.bam
	rm -f mapped/temp_${name}_merged.bam
	samtools index -@ $threads mapped/${name}_merged.bam
fi

### Mapping and identifying sRNA loci with shortstack
if [ ! -s mapped/${name}/ShortStack_All.gff3 ]; then
	rm -fr mapped/${name}
	ShortStack --bamfile mapped/${name}_merged.bam --genomefile ${ref_dir}/$ref.fa --sort_mem 8G --mmap u --dicermin 20 --dicermax 24 --bowtie_m all --mismatches 1 --foldsize 1000 --pad 250 --outdir mapped/${name}
fi

### Making bigiwig tracks
if [ ! -s tracks/${name}_merged_plus.bw ]; then
	printf "\nMaking tracks for $name\n"
	bamCoverage --filterRNAstrand forward -bs 1 -p $threads --normalizeUsing CPM -b mapped/${name}_merged.bam -o tracks/${name}_merged_plus.bw
	bamCoverage --filterRNAstrand reverse -bs 1 -p $threads --normalizeUsing CPM -b mapped/${name}_merged.bam -o tracks/${name}_merged_minus.bw
fi

touch chkpts/${name}
printf "\nsRNA mapping script finished successfully for $name\n"
