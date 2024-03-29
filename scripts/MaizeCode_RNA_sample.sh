#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=8G
#$ -l tmp_free=16G
#$ -o RNAsample.log
#$ -j y
#$ -N RNAsample

usage="
##### Script for Maize code RNA data analysis, used by script MaizeCode.sh with RNA argument
#####
##### sh MaizeCode_RNA_sample.sh -x data -d reference directory -l inbred line -t tissue -m RNA -r replicate ID -i sample ID -f path to sample -p paired -s step -a mappingoption
##### 	-x: type of data [ RNAseq | RAMPAGE ]
##### 	-d: folder containing the reference directory (e.g. ~/data/Genomes/Zea_mays/B73_v4)
##### 	-l: sample line (e.g. B73)
##### 	-t: tissue (e.g. endosperm)
##### 	-m: type of RNA sample [ RNAseq | RAMPAGE ]
##### 	-r: replicate ID (e.g. Rep1)
#####	-i: sample ID (name in original folder or SRR number)
#####	-f: path to original folder or SRA
##### 	-p: if data is paired-end (PE) or single-end (SE)
#####	-s: status of the raw data [ download | trim | done ] 'download' if sample needs to be copied/downloaded, 'trim' if only trimming has to be performed, 'done' if trimming has already been performed
#####	-a: what option to use for mapping [ default | all | colcen | colcenall ] (colcen: very-sensitive, -k 100; all: no MAPQ>10)
##### 	-h: help, returns usage
#####
##### It downloads or copies the files, runs fastQC, trims adapters with cutadapt, aligns with STAR (different parameters based on type of RNA),
##### make stranded track files (bigwigs) and get some mapping stats
#####
##### Requirements: samtools, fastQC, Cutadapt, STAR, Deeptools, parallel-fastq-dump (if downloading SRA data)
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
		m)	export rnatype=${OPTARG};;
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

if [ ! ${data} ] || [ ! ${ref_dir} ] || [ ! ${line} ] || [ ! ${tissue} ] || [ ! ${rnatype} ] || [ ! ${rep} ] || [ ! ${sampleID} ] || [ ! ${path} ] || [ ! ${paired} ] || [ ! ${step} ]; then
	printf "Missing arguments!\n"
	printf "${usage}\n"
	exit 1
fi

if [ ! ${mapparam} ]; then
	printf "No mapping option selected, using default\n"
	export mapparam="default"
elif [[ "${mapparam}" == "default" || "${mapparam}" == "colcen" || "${mapparam}" == "colcenall" || "${mapparam}" == "all" ]]; then
	printf "${mapparam} chosen as the mapping option\n"
else
	printf "Unknown mapping option selected\n"
	printf "${usage}\n"
	exit 1
fi

export ref=${ref_dir##*/}

name=${line}_${tissue}_${rnatype}_${rep}

case "${rnatype}" in
	RNAseq)		param_bg="--outWigType bedGraph"
			filesorder="fastq/trimmed_${name}_R1.fastq.gz fastq/trimmed_${name}_R2.fastq.gz"
			strandedness="reverse";;
	RAMPAGE) 	param_bg="--outWigType bedGraph read1_5p"
			filesorder="fastq/trimmed_${name}_R2.fastq.gz fastq/trimmed_${name}_R1.fastq.gz"
			strandedness="forward";;				
esac
	
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
		cutadapt -j ${threads} -q 10 -m 15 -a AGATCGGAAGAGCACACGTCTGAAC -A AGATCGGAAGAGCGTCGTGTAGGGA -o fastq/trimmed_${name}_R1.fastq.gz -p fastq/trimmed_${name}_R2.fastq.gz fastq/${name}_R1.fastq.gz fastq/${name}_R2.fastq.gz |& tee reports/trimming_${name}.txt
		#### Removing untrimmed fastq
		rm -f fastq/${name}_R*.fastq.gz
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for ${name}\n"
		fastqc -o reports/ fastq/trimmed_${name}_R1.fastq.gz
		fastqc -o reports/ fastq/trimmed_${name}_R2.fastq.gz
	fi
	#### Aligning reads to reference genome with STAR
	printf "\nMapping ${name} to ${ref} with STAR version:\n"
	STAR --version
	STAR --runMode alignReads --genomeDir ${ref_dir}/STAR_index --readFilesIn ${filesorder} --readFilesCommand zcat --runThreadN ${threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix mapped/map_${name}_ --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 20 --quantMode GeneCounts	
	### Marking duplicates
	STAR --runMode inputAlignmentsFromBAM --inputBAMfile mapped/map_${name}_Aligned.sortedByCoord.out.bam --bamRemoveDuplicatesType UniqueIdentical --outFileNamePrefix mapped/mrkdup_${name}_
	#### Indexing bam file
	printf "\nIndexing bam file\n"
	samtools index -@ ${threads} mapped/mrkdup_${name}_Processed.out.bam
	#### Getting stats from bam file
	printf "\nGetting some stats\n"
	samtools flagstat -@ ${threads} mapped/mrkdup_${name}_Processed.out.bam > reports/flagstat_${name}.txt
	### Making BedGraph files
	printf "\nMaking bedGraph files\n"
	STAR --runMode inputAlignmentsFromBAM --inputBAMfile mapped/mrkdup_${name}_Processed.out.bam --outWigStrand Stranded ${param_bg} --outFileNamePrefix tracks/bg_${name}_
	### Converting to bigwig files
	printf "\nConverting bedGraphs to bigWigs\n"
	bedSort tracks/bg_${name}_Signal.UniqueMultiple.str1.out.bg tracks/${name}_Signal.sorted.UniqueMultiple.str1.out.bg
	bedSort tracks/bg_${name}_Signal.Unique.str1.out.bg tracks/${name}_Signal.sorted.Unique.str1.out.bg
	bedSort tracks/bg_${name}_Signal.UniqueMultiple.str2.out.bg tracks/${name}_Signal.sorted.UniqueMultiple.str2.out.bg
	bedSort tracks/bg_${name}_Signal.Unique.str2.out.bg tracks/${name}_Signal.sorted.Unique.str2.out.bg
	if [[ ${strandedness} == "forward" ]]; then
		bedGraphToBigWig tracks/${name}_Signal.sorted.UniqueMultiple.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_plus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.Unique.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_unique_plus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.UniqueMultiple.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_minus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.Unique.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_unique_minus.bw
	elif [[ ${strandedness} == "reverse" ]]; then
		bedGraphToBigWig tracks/${name}_Signal.sorted.UniqueMultiple.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_minus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.Unique.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_unique_minus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.UniqueMultiple.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_plus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.Unique.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_unique_plus.bw
	else
		printf "\nStrandedness of data unknown! Tracks could not be created\n"
		exit 1
	fi	
	### Moving Log files to report folder
	mv mapped/*${name}*Log* reports/
	mv tracks/*${name}*Log* reports/
	### Cleaning up
	rm -f tracks/*${name}_Signal*
	### Summary stats
	printf "\nMaking mapping statistics summary\n"
	tot=$(grep "Total read pairs processed:" reports/trimming_${name}.txt | awk '{print $NF}' | sed 's/,//g')
	filt=$(grep "Number of input reads" reports/map_${name}_Log.final.out | awk '{print $NF}')
	multi=$(grep "Number of reads mapped to multiple loci" reports/map_${name}_Log.final.out | awk '{print $NF}')
	single=$(grep "Uniquely mapped reads number" reports/map_${name}_Log.final.out | awk '{print $NF}')
	allmap=$((multi+single))
	awk -v OFS="\t" -v l=${line} -v t=${tissue} -v m=${rnatype} -v r=${rep} -v g=${ref} -v a=${tot} -v b=${filt} -v c=${allmap} -v d=${single} 'BEGIN {print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}' >> reports/summary_mapping_stats.txt
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
		cutadapt -j ${threads} -q 10 -m 15 -a AGATCGGAAGAGCACACGTCTGAAC -o fastq/trimmed_${name}.fastq.gz fastq/${name}.fastq.gz |& tee reports/trimming_${name}.txt
		#### Removing untrimmed fastq
		rm -f fastq/${name}.fastq.gz
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for ${name}\n"
		fastqc -o reports/ fastq/trimmed_${name}.fastq.gz
	fi
	#### Aligning reads to reference genome with STAR
	printf "\nMapping ${name} to ${ref} with STAR version:\n"
	STAR --version
	STAR --runMode alignReads --genomeDir ${ref_dir}/STAR_index --readFilesIn fastq/trimmed_${name}.fastq.gz --readFilesCommand zcat --runThreadN ${threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix mapped/map_${name}_ --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 ${param_map} --quantMode GeneCounts |& tee reports/mapping_${name}.txt	
	#### Indexing bam file
	printf "\nIndexing bam file\n"
	samtools index -@ ${threads} mapped/map_${name}_Aligned.sortedByCoord.out.bam
	#### Getting stats from bam file
	printf "\nGetting some stats\n"
	samtools flagstat -@ ${threads} mapped/map_${name}_Aligned.sortedByCoord.out.bam > reports/flagstat_${name}.txt
	### Making BedGraph files
	printf "\nMaking bedGraph files\n"
	STAR --runMode inputAlignmentsFromBAM --inputBAMfile mapped/map_${name}_Aligned.sortedByCoord.out.bam --outWigStrand Stranded ${param_bg} --outFileNamePrefix tracks/bg_${name}_
	### Converting to bigwig files	
	printf "\nConverting bedGraphs to bigWigs\n"
	bedSort tracks/bg_${name}_Signal.UniqueMultiple.str1.out.bg tracks/${name}_Signal.sorted.UniqueMultiple.str1.out.bg
	bedSort tracks/bg_${name}_Signal.Unique.str1.out.bg tracks/${name}_Signal.sorted.Unique.str1.out.bg
	bedSort tracks/bg_${name}_Signal.UniqueMultiple.str2.out.bg tracks/${name}_Signal.sorted.UniqueMultiple.str2.out.bg
	bedSort tracks/bg_${name}_Signal.Unique.str2.out.bg tracks/${name}_Signal.sorted.Unique.str2.out.bg
	if [[ ${strandedness} == "forward" ]]; then
		bedGraphToBigWig tracks/${name}_Signal.sorted.UniqueMultiple.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_plus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.Unique.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_unique_plus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.UniqueMultiple.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_minus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.Unique.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_unique_minus.bw
	elif [[ ${strandedness} == "reverse" ]]; then
		bedGraphToBigWig tracks/${name}_Signal.sorted.UniqueMultiple.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_minus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.Unique.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_unique_minus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.UniqueMultiple.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_plus.bw
		bedGraphToBigWig tracks/${name}_Signal.sorted.Unique.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_unique_plus.bw
	else
		printf "\nStrandedness of data unknown! Tracks could not be created\n"
		exit 1
	fi	
	### Moving Log files to report folder
	mv mapped/*${name}*Log* reports/
	mv tracks/*${name}*Log* reports/
	### Cleaning up
	rm -f tracks/*${name}_Signal*
	### Summary stats
	printf "\nMaking mapping statistics summary\n"
	tot=$(grep "Total reads processed:" reports/trimming_${name}.txt | awk '{print $NF}' | sed 's/,//g')
	filt=$(grep "Number of input reads" reports/map_${name}_Log.final.out | awk '{print $NF}')
	multi=$(grep "Number of reads mapped to multiple loci" reports/map_${name}_Log.final.out | awk '{print $NF}')
	single=$(grep "Uniquely mapped reads number" reports/map_${name}_Log.final.out | awk '{print $NF}')
	allmap=$((multi+single))
	awk -v OFS="\t" -v l=${line} -v t=${tissue} -v m=${rnatype} -v r=${rep} -v g=${ref} -v a=${tot} -v b=${filt} -v c=${allmap} -v d=${single} 'BEGIN {print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}' >> reports/summary_mapping_stats.txt
else
	printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
	exit 1
fi


printf "\nScript finished successfully!\n"
touch chkpts/${name}_${ref}

