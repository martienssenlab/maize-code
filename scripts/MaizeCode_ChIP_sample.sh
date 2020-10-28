#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=12G
#$ -l tmp_free=100G
#$ -o logs/ChIPsample.log
#$ -j y
#$ -N ChIPsample

usage="
##### Script for Maize code ChIP data analysis, used by script MaizeCode.sh for ChIP samples
#####
##### sh MaizeCode_ChIP_sample.sh -d reference directory -l inbred line -t tissue -m histone mark -e replicate ID -p paired
##### 	-d: folder containing the reference directory (e.g. ~/data/Genomes/Zea_mays/B73_v4)
##### 	-l: inbred line (e.g. B73)
##### 	-t: tissue (e.g. endosperm)
##### 	-m: ChIP-seq mark (e.g. H3K4me1)
##### 	-e: replicate ID (e.g. Rep1)
##### 	-p: if data is paired-end (PE) or single-end (SE) [ PE | SE ]
##### 	-h: help, returns usage
#####
##### It runs fastQC, trims adapters with cutadapt, aligns with bowtie2,
##### filters duplicates with samtools, and get some mapping stats
##### !!! At the moment, it is not going to be run on the same references, but will overwrite if the same sample name is mapping to a new reference !!!
#####
##### Requirements: samtools, fastQC, Cutadapt, Bowtie2
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

while getopts "d:l:t:m:r:p:h" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		d) 	export ref_dir=${OPTARG};;
		l)	export line=${OPTARG};;
		t)	export tissue=${OPTARG};;
		m)	export mark=${OPTARG};;
		r)	export rep=${OPTARG};;
		p)	export paired=${OPTARG};;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $ref_dir ] || [ ! $line ] || [ ! $tissue ] || [ ! $mark ] || [ ! $rep ] || [ ! $paired ]; then
	printf "Missing arguments!\n"
	printf "$usage\n"
	exit 1
fi

export ref=${ref_dir##*/}

name=${line}_${tissue}_${mark}_${rep}

if [[ $paired == "PE" ]]; then
	#### FastQC on raw data
	printf "\nRunning fastQC for $name with fastqc version:\n"
	fastqc --version
	fastqc -o reports/ fastq/${name}_R1.fastq.gz
	fastqc -o reports/ fastq/${name}_R2.fastq.gz	
	#### Trimming illumina adapters with Cutadapt
	printf "\nTrimming Illumina adapters for $name with cutadapt version:\n"
	cutadapt --version
	cutadapt -j $threads -q 10 -m 20 -a AGATCGGAAGAGCACACGTCTGAAC -A AGATCGGAAGAGCGTCGTGTAGGGA -o fastq/trimmed_${name}_R1.fastq.gz -p fastq/trimmed_${name}_R2.fastq.gz fastq/${name}_R1.fastq.gz fastq/${name}_R2.fastq.gz |& tee reports/trimming_${name}.txt
	#### FastQC on trimmed data
	printf "\nRunning fastQC on trimmed files for $name\n"
	fastqc -o reports/ fastq/trimmed_${name}_R1.fastq.gz
	fastqc -o reports/ fastq/trimmed_${name}_R2.fastq.gz
	#### Aligning reads to reference genome with Bowtie2
	#### maxins 1500 used after seeing that average insert size from first round of mapping was ~500bp (for most B73 marks) but ~900bp for Inputs
	printf "\nMaping $name to $ref\n"
	bowtie2 --version
	bowtie2 -p $threads --end-to-end --maxins 1500 --met-file reports/bt2_${name}.txt -x $ref_dir/$ref -1 fastq/trimmed_${name}_R1.fastq.gz -2 fastq/trimmed_${name}_R2.fastq.gz -S mapped/${name}.sam |& tee reports/mapping_${name}.txt
elif [[ $paired == "SE" ]]; then
### Single-end process not tested !!!
	#### FastQC on raw data
	printf "\nRunning fastQC for $name with fastqc version:\n"
	fastqc --version
	fastqc -o reports/ fastq/${name}.fastq.gz
	#### Trimming illumina adapters with Cutadapt
	printf "\nTrimming Illumina adapters for $name with cutadapt version:\n"
	cutadapt --version
	cutadapt -j $threads -q 10 -m 20 -a AGATCGGAAGAGCACACGTCTGAAC -o fastq/trimmed_${name}.fastq.gz fastq/${name}.fastq.gz |& tee reports/trimming_${name}.txt
	#### FastQC on trimmed data
	printf "\nRunning fastQC on trimmed files for $name\n"
	fastqc -o reports/ fastq/trimmed_${name}.fastq.gz
	#### Aligning reads to reference genome with Bowtie2
	printf "\nMaping $name to $ref with bowtie2 version:\n"
	bowtie2 --version
	bowtie2 -p $threads --end-to-end --met-file reports/bt2_${name}.txt -x $ref_dir/$ref -U fastq/trimmed_${name}.fastq.gz -S mapped/${name}.sam |& tee reports/mapping_${name}.txt
else
	printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
	exit 1
fi

#### Removing duplicates, sorting, converting to bam and indexing file with samtools
printf "\nRemoving duplicates, sorting and indexing file with samtools version:\n"
samtools --version
samtools fixmate -@ $threads -m mapped/${name}.sam mapped/temp1_${name}.bam
samtools sort -@ $threads -o mapped/${name}.bam mapped/temp1_${name}.bam
rm -f mapped/${name}.sam
rm -f mapped/temp1_${name}.bam
samtools markdup -r -s -f reports/markdup_${name}.txt -@ $threads mapped/${name}.bam mapped/rmdup_${name}.bam
samtools index -@ $threads mapped/rmdup_${name}.bam
printf "\nGetting some stats\n"
samtools flagstat -@ $threads mapped/rmdup_${name}.bam > reports/flagstat_${name}.txt

#### Summary stats only working with PE for now
if [[ $paired == "PE" ]]; then
	printf "\nMaking mapping statistics summary\n"
	totp=$(grep "Total read pairs processed:" reports/trimming_${name}.txt | awk '{print $NF}' | sed 's/,//g')
	filtp=$(grep "reads" reports/mapping_${name}.txt | awk '{print $1}')
	dedup=$(grep "total" reports/flagstat_${name}.txt | awk '{print $1}')
	prop=$(grep "properly paired" reports/flagstat_${name}.txt | awk '{print $1}')
	tots=$((totp*2))
	filts=$((filtp*2))
	awk -v OFS="\t" -v l=$line -v t=$tissue -v m=$mark -v r=$rep -v a=$tots -v b=$filts -v c=$dedup -v d=$prop 'BEGIN {print l,t,m,r,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}' >> reports/summary_mapping_stats.txt
fi

printf "\nScript finished successfully!\n"
touch chkpts/${name}_${ref}
