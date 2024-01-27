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
#####	-s: status of the raw data [ download | trim | qc | align | filter | done ] 'download' if sample needs to be copied/downloaded, 'trim' if only trimming has to be performed, etc. 'done' if it has already been processed to aligned, filtered BAMs.
#####	-a: what option to use for mapping [ default | all | colcen | colcenall ] (colcen: very-sensitive, -k 100; all: no MAPQ>10)
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

if [ ! ${data} ] || [ ! ${ref_dir} ] || [ ! ${line} ] || [ ! ${tissue} ] || [ ! ${mark} ] || [ ! ${rep} ] || [ ! ${sampleID} ] || [ ! ${path} ] || [ ! ${paired} ] || [ ! ${step} ]; then
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
	
if [[ ! (${paired} == "PE" || ${paired} == "SE") ]]; then
	printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
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


if [[ ${step} == "download" ]]; then
	if [[ ${path} == "SRA" ]]; then
		printf "\nUsing fasterq-dump for ${name} (${sampleID})\n"
		fasterq-dump -e ${threads} --outdir ./fastq ${sampleID}
		printf "\n$name ($sampleID) downloaded\nGzipping and renaming files..."

		if [[ ${paired} == "PE" ]]; then
			pigz -p ${threads} ./fastq/${sampleID}_1.fastq
			mv ./fastq/${sampleID}_1.fastq.gz ./fastq/${name}_R1.fastq.gz
			pigz -p ${threads} ./fastq/${sampleID}_2.fastq
			mv ./fastq/${sampleID}_2.fastq.gz ./fastq/${name}_R2.fastq.gz
		elif [[ ${paired} == "SE" ]]; then
			pigz -p ${threads} ./fastq/${sampleID}.fastq
			mv ./fastq/${sampleID}.fastq.gz ./fastq/${name}.fastq.gz
		fi
	else
		printf "\nCopying ${paired} fastq(s) for ${name} (${sampleID} in ${path})\n"
		if [[ ${paired} == "PE" ]]; then
			cp ${path}/*${sampleID}*R1*q.gz ./fastq/${name}_R1.fastq.gz
			cp ${path}/*${sampleID}*R2*q.gz ./fastq/${name}_R2.fastq.gz
		elif [[ ${paired} == "SE" ]]; then
			cp ${path}/*${sampleID}*q.gz ./fastq/${name}.fastq.gz
		fi
	fi

	step="trim"
fi


#### Trimming illumina adapters with Cutadapt
if [[ ${step} == "trim" ]]; then
	printf "\nTrimming Illumina adapters for ${name} with cutadapt version:\n"
	cutadapt --version

	if [[ ${paired} == "PE" ]]; then
		ca_adapter_params="-a AGATCGGAAGAGCACACGTCTGAAC -A AGATCGGAAGAGCGTCGTGTAGGGA"
		ca_read_params="-o fastq/trimmed_${name}_R1.fastq.gz -p fastq/trimmed_${name}_R2.fastq.gz fastq/${name}_R1.fastq.gz fastq/${name}_R2.fastq.gz"

	elif [[ ${paired} == "SE" ]]; then
		ca_adapter_params="-a AGATCGGAAGAGCACACGTCTGAAC"
		ca_read_params="-o fastq/trimmed_${name}.fastq.gz fastq/${name}.fastq.gz"
	fi
	
	cutadapt -j ${threads} -q 10 -m 20 ${ca_adapter_params} ${ca_read_params} |& tee reports/trimming_${name}.txt

	#### Removing untrimmed fastq(s)
	if [[ ${paired} == "PE" ]]; then
		rm -f fastq/${name}_R{1,2}.fastq.gz

	elif [[ ${paired} == "SE" ]]; then
		rm -f fastq/${name}.fastq.gz
	fi

	step="qc"
fi


#### FastQC on raw and trimmed reads
if [[ ${step} == "qc" ]]; then
	printf "\nRunning fastQC on raw and trimmed reads for ${name} with fastqc version:\n"
	fastqc --version
	if [[ ${paired} == "PE" ]]; then
		fastqc --threads ${threads} -o reports/ fastq/${name}_R1.fastq.gz fastq/${name}_R2.fastq.gz fastq/trimmed_${name}_R1.fastq.gz fastq/trimmed_${name}_R2.fastq.gz
	elif [[ ${paired} == "SE" ]]; then
		fastqc --threads ${threads} -o reports/ fastq/${name}.fastq.gz fastq/trimmed_${name}.fastq.gz
	fi

	step="align"
fi


if [[ ${step} == "align" ]]; then
	#### Align reads to the reference respecting the sequencing strategy and mapping macro parameter
	printf "\nMapping ${paired} library ${name} to ${ref} with ${mapparam} parameters with bowtie2 version:\n"
	bowtie2 --version
	printf "\nSecondary alignments will be discarded.\n"

	if [[ ${paired} == "PE" ]]; then
		#### maxins 1500 used after seeing that average insert size from first round of mapping was ~500bp (for most B73 marks) but ~900bp for Inputs
		paired_params="--maxins 1500 -1 fastq/trimmed_${name}_R1.fastq.gz -2 fastq/trimmed_${name}_R2.fastq.gz"
	elif [[ ${paired} == "SE" ]]; then
		paired_params="-U fastq/trimmed_${name}.fastq.gz"
	fi

	if [[ ${mapparam} == "default" || ${mapparam} == "all" ]]; then
		map_params=""
	elif [[ ${mapparam} == "colcen" || ${mapparam} == "colcenall" ]]; then
		map_params="--very-sensitive --no-mixed --no-discordant --k 100"
	fi
	
	step="filter"
fi

# The shell redirection to a third file descriptor here is necessary only because we tee the bt2 stderr info into the log and the mapping report.
# Otherwise, could just write the report directly to the with with 2>reports/mapping_${name}.txt and pipe stdout to samtools.
( bowtie2 -p ${threads} --end-to-end --met-file reports/bt2_${name}.txt -x $ref_dir/$ref ${paired_params} ${map_params} \
	| samtools view -@ ${threads} -b -h -F 256 -o mapped/temp1_${name}.bam) 3>&1 1>&2 2>&3 | tee reports/mapping_${name}.txt
	

if [[ ${step} == "filter" ]]; then
	#### Filter read alignments respecting the sequencing strategy and mapping macro parameter
	printf "\nRemoving low quality reads, secondary alignements and duplicates, sorting and indexing file with samtools version:\n"
	samtools --version

	if [[ ${mapparam} == "default" || ${mapparam} == "colcen" ]]; then
		filter_params="-q10"
	elif [[ ${mapparam} == "colcenall" || ${mapparam} == "all" ]]; then
		filter_params=""
	fi

	samtools view -@ ${threads} -b -h -u ${filter_params} mapped/temp1_${name}.bam \
		| samtools fixmate -@ ${threads} -m -u - - \
		| samtools sort -@ ${threads} -m 1G -u -T temp1_${name} - \
		| samtools markdup -r -s -f reports/markdup_${name}.txt -@ ${threads} - mapped/${name}.bam
	samtools index -@ ${threads} mapped/${name}.bam
fi


#### Cleanup
rm -f mapped/temp*_${name}.bam


#### Summary stats
printf "\nGetting some stats\n"
samtools flagstat -@ ${threads} mapped/${name}.bam > reports/flagstat_${name}.txt

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
awk -v OFS="\t" \
	-v l=${line} \
	-v t=${tissue} \
	-v m=${mark} \
	-v r=${rep}${add} \
	-v g=${ref} \
	-v a=${tot} \
	-v b=${filt} \
	-v c=${allmap} \
	-v d=${single} \
	'BEGIN {print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}' \
	>> reports/summary_mapping_stats.txt

printf "\nScript finished successfully!\n"
touch chkpts/${name}_${ref}

