#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=6G
#$ -l tmp_free=8G
#$ -o logs/shRNAsample.log
#$ -j y
#$ -N shRNAsample

usage="
##### Script for Maize code shRNA data analysis, used by script MaizeCode.sh for shRNA samples
#####
##### sh MaizeCode_shRNA_sample.sh -x datatype -d reference directory -l inbred line -t tissue -m shRNA -r replicate ID -i sample ID -f path to sample -p paired -s step
##### 	-x: type of data (not used here yet mandatory, should be 'shRNA')
##### 	-d: folder containing the reference directory (e.g. ~/data/Genomes/Zea_mays/B73_v4)
##### 	-l: inbred line (e.g. B73)
##### 	-t: tissue (e.g. endosperm)
##### 	-m: sample name (should only be 'shRNA' for now at least)
##### 	-r: replicate ID (e.g. Rep1)
#####	  -i: sample ID (name in original folder or SRR number)
#####	  -f: path to original folder or SRA
##### 	-p: if data is paired-end (PE) or single-end (SE) [ PE | SE ]
#####	  -s: status of the raw data [ download | trim | done ] 'download' if sample needs to be copied/downloaded, 'trim' if only trimming has to be performed, 'done' if trimming has already been performed
##### 	-h: help, returns usage
#####
##### It downloads or copies the files, runs fastQC, trims adapters with cutadapt, aligns to structural RNA with bowtie2 to filter them out,
##### aligns the rest and call sRNA clusters with Shortstack, creates bigiwig files with deeptools, filters sizes with seqkit and get some mapping stats
#####
##### Requirements: samtools, fastQC, Cutadapt, Bowtie2, ShortStack, seqkit, parallel-fastq-dump (if downloading SRA data)
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=${NSLOTS}

if [ $# -eq 0 ]; then
	printf "${usage}\n"
	exit 1
fi

while getopts "x:d:l:t:m:r:i:f:p:s:h" opt; do
	case ${opt} in
		h) 	printf "${usage}\n"
			exit 0;;
		x)	export data=${OPTARG};;
		d) 	export ref_dir=${OPTARG};;
		l)	export line=${OPTARG};;
		t)	export tissue=${OPTARG};;
		m)	export rnatype=${OPTARG};;
		r)	export rep=${OPTARG};;
		i)	export sampleID=${OPTARG};;
		f)	export path=${OPTARG};;
		p)	export paired=${OPTARG};;
		s)	export step=${OPTARG};;
		*)	printf "${usage}\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! ${data} ] || [ ! ${ref_dir} ] || [ ! ${line} ] || [ ! ${tissue} ] || [ ! ${rnatype} ] || [ ! ${rep} ] || [ ! ${sampleID} ] || [ ! ${path} ] || [ ! ${paired} ] || [ ! ${step} ]; then
	printf "Missing arguments!\n"
	printf "${usage}\n"
	exit 1
fi

export ref=${ref_dir##*/}

name=${line}_${tissue}_${rnatype}_${rep}

if [ -s ${ref_dir}/*.fa ]; then
	fa_file=$(ls ${ref_dir}/*.fa)
	fasta=${fa_file}
elif [ -s ${ref_dir}/*.fasta ]; then
	fa_file=$(ls ${ref_dir}/*.fasta)
	fasta=${fa_file}
elif [ -s ${ref_dir}/*.fa.gz ]; then
	fa_file=$(ls ${ref_dir}/*.fa.gz)
	pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${data}_${ref}.fa
	fasta=${ref_dir}/temp_${data}_${ref}.fa
elif [ -s ${ref_dir}/*.fasta.gz ]; then
	fa_file=$(ls ${ref_dir}/*.fasta.gz)
	pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${data}_${ref}.fa
	fasta=${ref_dir}/temp_${data}_${ref}.fa
fi

if [[ ${paired} == "PE" ]]; then
	
	###### THIS PART OF THE PE PIPELINE WILL NOT WORK AS IS !!! NOT THAT USEFUL ANYWAY, BUT NEED TO BE CHANGED JUST IN CASE ######
	###### Need to change adapter sequences based on kit used !! And how to process PE files after filtering out with bowtie2 ####
	
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
		cutadapt -j ${threads} -q 10 -m 15 -a AGATCGGAAGAGCACACGTCTGAAC -A AGATCGGAAGAGCGTCGTGTAGGGA -o fastq/trimmed_${name}_R1.fastq.gz -p fastq/trimmed_${name}_R2.fastq.gz fastq/${name}_R1.fastq.gz fastq/${name}_R2.fastq.gz |& tee reports/trimming_${name}.txt
		#### Removing untrimmed fastq
		rm -f fastq/${name}_R*.fastq.gz
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for ${name}\n"
		fastqc -o reports/ fastq/trimmed_${name}_R1.fastq.gz
		fastqc -o reports/ fastq/trimmed_${name}_R2.fastq.gz
	fi
  	#### Aligning reads to filter out structural RNAs (rRNAs, snoRNAs and tRNAs) with bowtie2
	bowtie2 --very-sensitive -p ${threads} -x structural_RNA/zm_structural_RNAs -1 fastq/trimmed_${name}_R1.fastq.gz -2 fastq/trimmed_${name}_R2.fastq.gz | samtools view -@ ${threads} -f 0x4 | samtools fastq -@ ${threads} | gzip > fastq/filtered_${name}.fastq.gz
	#### Mapping and identifying sRNA loci with shortstack
	ShortStack --readfile fastq/filtered_${name}.fastq.gz --genomefile ${fasta} --bowtie_cores $threads --sort_mem 4G --mmap u --dicermin 20 --dicermax 24 --bowtie_m all --mismatches 1 --foldsize 1000 --pad 250 --outdir mapped/${name}
	#### Making bigiwig tracks
	samtools index -@ $threads mapped/${name}/filtered_${name}.bam
	printf "\nMaking plus track for $name\n"
	bamCoverage --filterRNAstrand forward -bs 1 -p $threads --normalizeUsing CPM -b mapped/${name}/filtered_${name}.bam -o tracks/${name}_plus.bw
	bamCoverage --filterRNAstrand reverse -bs 1 -p $threads --normalizeUsing CPM -b mapped/${name}/filtered_${name}.bam -o tracks/${name}_minus.bw
	#### Filtering only small RNA sizes (15 to 32nt)
#	seqkit seq --max-len 32 --min-len 15 fastq/filtered_${name}.fastq.gz | gzip > fastq/sized_${name}.fastq.gz
	touch chkpts/${name}
	
	###########
  
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
		cutadapt -j ${threads} -q 10 -m 15 -a TGGAATTCTCGGGTGCCAAGG -o fastq/trimmed_${name}.fastq.gz fastq/${name}.fastq.gz |& tee reports/trimming_${name}.txt
		#### Removing untrimmed fastq
		rm -f fastq/${name}.fastq.gz
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for ${name}\n"
		fastqc -o reports/ fastq/trimmed_${name}.fastq.gz
	fi
	#### Aligning reads to filter out structural RNAs (rRNAs, snoRNAs and tRNAs) with bowtie2
  	printf "\nAligning reads to filter out structural RNAs with bowtie2 version:\n"
  	bowtie2 --version
	bowtie2 --very-sensitive -p ${threads} -x structural_RNA/zm_structural_RNAs -U fastq/trimmed_${name}.fastq.gz | samtools view -@ ${threads} -f 0x4 | samtools fastq -@ ${threads} | gzip > fastq/filtered_${name}.fastq.gz
	#### Mapping and identifying sRNA loci with shortstack
	ShortStack --readfile fastq/filtered_${name}.fastq.gz --genomefile ${fasta} --bowtie_cores ${threads} --sort_mem 6G --mmap u --dicermin 20 --dicermax 24 --bowtie_m 1000 --mismatches 1 --foldsize 1000 --pad 250 --outdir mapped/${name}
	#### Making bigiwig tracks
	samtools index -@ ${threads} mapped/${name}/filtered_${name}.bam
	printf "\nMaking plus track for ${name}\n"
	bamCoverage --filterRNAstrand forward -bs 1 -p ${threads} --normalizeUsing CPM -b mapped/${name}/filtered_${name}.bam -o tracks/${name}_plus.bw
	bamCoverage --filterRNAstrand reverse -bs 1 -p ${threads} --normalizeUsing CPM -b mapped/${name}/filtered_${name}.bam -o tracks/${name}_minus.bw
	#### Filtering only small RNA sizes (15 to 32nt)
	samtools view -h mapped/${name}/filtered_${name}.bam | awk '(length($10) >= 20 && length($10) <= 24) || $1 ~ /^@/' | samtools view -bS - > mapped/${name}/sized_${name}.bam
	#### Getting stats of size distribution
  	printf "\nGetting trimmed stats for ${name}\n"
  	zcat fastq/trimmed_${name}.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c | awk -v OFS="\t" -v n=${name} '{print n,"trimmed",$2,$1}' > reports/sizes_trimmed_${name}.txt
  	printf "\nGetting filtered stats for ${name}\n"
	zcat fastq/filtered_${name}.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c | awk -v OFS="\t" -v n=${name} '{print n,"filtered",$2,$1}' > reports/sizes_filtered_${name}.txt
  	printf "\nGetting mapped stats for ${name}\n"
	samtools view mapped/${name}/filtered_${name}.bam | awk '$2==0 || $2==16 {print length($10)}' | sort -n | uniq -c | awk -v OFS="\t" -v n=${name} '{print n,"mapped",$2,$1}' > reports/sizes_mapped_${name}.txt
else
	printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
	exit 1
fi

printf "\nScript finished successfully!\n"
touch chkpts/${name}_${ref}
