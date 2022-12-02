#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=1.5G
#$ -l tmp_free=20G
#$ -o logs/ChIPsample.log
#$ -j y
#$ -N ChIPsample

usage="
##### Script for Maize code Histone ChIP data analysis, used by script MaizeCode.sh for ChIP samples
#####
##### sh MaizeCode_ChIP_sample.sh -x datatype -d reference directory -l inbred line -t tissue -m histone mark -r replicate ID -i sample ID -f path to sample -p paired -s step -a mappingoption
##### 	-x: type of data (should be 'ChIP' or 'ChIP_*' where * is a index marking which input to use when several inputs are present)
##### 	-d: folder containing the reference directory (e.g. ~/data/Genomes/Zea_mays/B73_v4)
##### 	-l: inbred line (e.g. B73)
##### 	-t: tissue (e.g. endosperm)
##### 	-m: ChIP-seq mark (e.g. H3K4me1)
##### 	-r: replicate ID (e.g. Rep1)
#####	-i: sample ID (name in original folder or SRR number)
#####	-f: path to original folder or SRA
##### 	-p: if data is paired-end (PE) or single-end (SE) [ PE | SE ]
#####	-s: status of the raw data [ download | trim | done ] 'download' if sample needs to be copied/downloaded, 'trim' if only trimming has to be performed, 'done' if trimming has already been performed
#####	-a: what option to use for mapping [ default | Colcen ]
##### 	-h: help, returns usage
#####
##### It downloads or copies the files, runs fastQC, trims adapters with cutadapt, aligns with bowtie2,
##### filters duplicates with samtools, and get some mapping stats
#####
##### Requirements: samtools, fastQC, Cutadapt, Bowtie2, parallel-fastq-dump (if downloading SRA data)
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

while getopts "x:d:l:t:m:r:i:f:p:s:a:h" opt; do
	case ${opt} in
		h) 	printf "${usage}\n"
			exit 0;;
		x)	export data=${OPTARG};;
		d) 	export ref_dir=${OPTARG};;		
		l)	export line=${OPTARG};;
		t)	export tissue=${OPTARG};;
		m)	export mark=${OPTARG};;
		r)	export rep=${OPTARG};;
		i)	export sampleID=${OPTARG};;
		f)	export path=${OPTARG};;
		p)	export paired=${OPTARG};;
		s)	export step=${OPTARG};;
		a)	export mapparam=${OPTARG};;
		*)	printf "${usage}\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! ${data} ] || [ ! ${ref_dir} ] || [ ! ${line} ] || [ ! ${tissue} ] || [ ! ${mark} ] || [ ! ${rep} ] || [ ! ${sampleID} ] || [ ! ${path} ] || [ ! ${paired} ] || [ ! ${step} ] || [ ! ${mapparam} ]; then
	printf "Missing arguments!\n"
	printf "${usage}\n"
	exit 1
fi

export ref=${ref_dir##*/}

if [[ ${data} == "ChIP_"* ]] && [[ ${mark} == "Input" ]]; then
	tmp=${data#ChIP_}
	add="_${tmp}"
	name=${line}_${tissue}_${mark}_${rep}${add}
else
	name=${line}_${tissue}_${mark}_${rep}
fi

if [[ ${paired} == "PE" ]]; then
	if [[ ${step} == "download" ]]; then
		if [[ ${path} == "SRA" ]]; then
			printf "\nUsing fasterq-dump for ${name} (${sampleID})\n"
			fasterq-dump -e ${threads} --outdir ./fastq ${sampleID}
			printf "\n$name ($sampleID) downloaded\nGzipping and renaming files..."
			pigz -p ${threads} ./fastq/${sampleID}_1.fastq
			mv ./fastq/${sampleID}_1.fastq.gz ./fastq/${name}_R1.fastq.gz
			pigz -p ${threads} ./fastq/${sampleID}_2.fastq
			mv ./fastq/${sampleID}_2.fastq.gz ./fastq/${name}_R2.fastq.gz
			step="trim"
		else
			printf "\nCopying PE fastq for ${name} (${sampleID} in ${path})\n"
			cp ${path}/*${sampleID}*R1*q.gz ./fastq/${name}_R1.fastq.gz
			cp ${path}/*${sampleID}*R2*q.gz ./fastq/${name}_R2.fastq.gz
			step="trim"
		fi
	fi
	if [[ ${step} == "trim" ]]; then
		#### FastQC on raw data
		printf "\nRunning fastQC for ${name} with fastqc version:\n"
		fastqc --version
		fastqc -o reports/ fastq/${name}_R1.fastq.gz
		fastqc -o reports/ fastq/${name}_R2.fastq.gz	
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for ${name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j ${threads} -q 10 -m 20 -a AGATCGGAAGAGCACACGTCTGAAC -A AGATCGGAAGAGCGTCGTGTAGGGA -o fastq/trimmed_${name}_R1.fastq.gz -p fastq/trimmed_${name}_R2.fastq.gz fastq/${name}_R1.fastq.gz fastq/${name}_R2.fastq.gz |& tee reports/trimming_${name}.txt
		#### Removing untrimmed fastq
		rm -f fastq/${name}_R*.fastq.gz
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for ${name}\n"
		fastqc -o reports/ fastq/trimmed_${name}_R1.fastq.gz
		fastqc -o reports/ fastq/trimmed_${name}_R2.fastq.gz
	fi
	#### Aligning reads to reference genome with Bowtie2
	#### maxins 1500 used after seeing that average insert size from first round of mapping was ~500bp (for most B73 marks) but ~900bp for Inputs
	if [[ ${mapparam} == "default" ]]; then
		printf "\nMaping ${name} to ${ref} with ${mapparam} parameters\n"
		bowtie2 --version
		bowtie2 -p ${threads} --end-to-end --maxins 1500 --met-file reports/bt2_${name}.txt -x $ref_dir/$ref -1 fastq/trimmed_${name}_R1.fastq.gz -2 fastq/trimmed_${name}_R2.fastq.gz -S mapped/${name}.sam |& tee reports/mapping_${name}.txt
	elif [[ ${mapparam} == "Colcen" ]]; then
		printf "\nMaping ${name} to ${ref} with ${mapparam} parameters\n"
		bowtie2 --version
		bowtie2 -p ${threads} -k 150 --end-to-end --maxins 1500 --met-file reports/bt2_${name}.txt -x $ref_dir/$ref -1 fastq/trimmed_${name}_R1.fastq.gz -2 fastq/trimmed_${name}_R2.fastq.gz -S mapped/${name}.sam |& tee reports/mapping_${name}.txt
	fi
elif [[ ${paired} == "SE" ]]; then
	if [[ ${step} == "download" ]]; then
		if [[ ${path} == "SRA" ]]; then
			printf "\nUsing fasterq-dump for ${name} (${sampleID})\n"
			fasterq-dump -e ${threads} --outdir ./fastq ${sampleID}
			printf "\n$name ($sampleID) downloaded\nRenaming files..."
			pigz -p ${threads} ./fastq/${sampleID}.fastq
			mv ./fastq/${sampleID}.fastq.gz ./fastq/${name}.fastq.gz
			step="trim"
		else
			printf "\nCopying SE fastq for ${name} (${sampleID} in ${path})\n"
			cp ${path}/*${sampleID}*q.gz ./fastq/${name}.fastq.gz
			step="trim"
		fi
	fi
	if [[ ${step} == "trim" ]]; then
		#### FastQC on raw data
		printf "\nRunning fastQC for ${name} with fastqc version:\n"
		fastqc --version
		fastqc -o reports/ fastq/${name}.fastq.gz
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for ${name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j ${threads} -q 10 -m 20 -a AGATCGGAAGAGCACACGTCTGAAC -o fastq/trimmed_${name}.fastq.gz fastq/${name}.fastq.gz |& tee reports/trimming_${name}.txt
		#### Removing untrimmed fastq
		rm -f fastq/${name}.fastq.gz
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for ${name}\n"
		fastqc -o reports/ fastq/trimmed_${name}.fastq.gz
	fi
	#### Aligning reads to reference genome with Bowtie2
	if [[ ${mapparam} == "default" ]]; then
		printf "\nMaping ${name} to ${ref} with ${mapparam} parameters\n"
		bowtie2 --version
		bowtie2 -p ${threads} --end-to-end --met-file reports/bt2_${name}.txt -x $ref_dir/$ref -U fastq/trimmed_${name}.fastq.gz -S mapped/${name}.sam |& tee reports/mapping_${name}.txt
	elif [[ ${mapparam} == "Colcen" ]]; then
		printf "\nMaping ${name} to ${ref} with ${mapparam} parameters\n"
		bowtie2 --version
		bowtie2 -p ${threads} --very-sensitive --no-mixed --no-discordant --k 100 --end-to-end --met-file reports/bt2_${name}.txt -x $ref_dir/$ref -U fastq/trimmed_${name}.fastq.gz -S mapped/${name}.sam |& tee reports/mapping_${name}.txt
	fi
else
	printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
	exit 1
fi

#### Removing duplicates, sorting, converting to bam and indexing file with samtools
printf "\nRemoving duplicates, sorting and indexing file with samtools version:\n"
samtools --version
samtools fixmate -@ ${threads} -m mapped/${name}.sam mapped/temp1_${name}.bam
rm -f mapped/${name}.sam
samtools sort -@ ${threads} -o mapped/temp2_${name}.bam mapped/temp1_${name}.bam
samtools markdup -r -s -S -f reports/markdup_${name}.txt -@ ${threads} mapped/temp2_${name}.bam mapped/${name}.bam
samtools index -@ ${threads} mapped/${name}.bam
printf "\nGetting some stats\n"
samtools flagstat -@ ${threads} mapped/${name}.bam > reports/flagstat_${name}.txt
rm -f mapped/temp*_${name}.bam
#### Summary stats
printf "\nMaking mapping statistics summary\n"
if [[ ${paired} == "PE" ]]; then
	tot=$(grep "Total read pairs processed:" reports/trimming_${name}.txt | awk '{print $NF}' | sed 's/,//g')
	filt=$(grep "reads" reports/mapping_${name}.txt | awk '{print $1}')
	multi=$(grep "aligned concordantly >1 times" reports/mapping_${name}.txt | awk '{print $1}')
	single=$(grep "aligned concordantly exactly 1 time" reports/mapping_${name}.txt | awk '{print $1}')
else
	tot=$(grep "Total reads processed:" reports/trimming_${name}.txt | awk '{print $NF}' | sed 's/,//g')
	filt=$(grep "reads" reports/mapping_${name}.txt | awk '{print $1}')
	multi=$(grep "aligned >1 times" reports/mapping_${name}.txt | awk '{print $1}')
	single=$(grep "aligned exactly 1 time" reports/mapping_${name}.txt | awk '{print $1}')
fi
allmap=$((multi+single))
awk -v OFS="\t" -v l=${line} -v t=${tissue} -v m=${mark} -v r=${rep}${add} -v g=${ref} -v a=${tot} -v b=${filt} -v c=${allmap} -v d=${single} 'BEGIN {print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}' >> reports/summary_mapping_stats.txt

printf "\nScript finished successfully!\n"
touch chkpts/${name}_${ref}
