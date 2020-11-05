#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=10G
#$ -l tmp_free=50G
#$ -o env.log
#$ -j y
#$ -N env

usage="
##### Script for checking environment and building indexes for Maize code
##### 
##### sh MaizeCode_check_environment.sh -p path to genome reference -r ref -d datatype
##### 	-p: path to the folder containing all the different genome references (e.g. ~/data/Genomes/Zea_mays)
##### 	-r: genome reference to use (e.g. B73_v4) 
#####	-d: type of data [ChIP | RNA]
##### 	-h: help, returns usage
#####
##### The reference genome folder should contain a single fasta file and a single gff3 file (can be gzipped)
##### The gff3 files should have 'gene' in column 3 and exons should be linked by 'Parent' in column 9
##### The fasta and gff3 files should have the same chromosome names (i.e. 1 2 3... and 1 2 3... or Chr1 Chr2 Chr3... and Chr1 Chr2 Chr3...)
##### The type of data (ChIP or RNA) will define whether to build bowtie2 or STAR indexes
#####
##### This script check if the indexes and other files need to be created and creates them if they do
#####
##### Requirements: pigz, samtools, Bowtie2 for ChIP data, STAR for RNA data
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

while getopts "p:r:d:h" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		p)	export pathtoref=${OPTARG};;
		r)	export ref=${OPTARG};;
		d)	export datatype=${OPTARG};;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $ref ] || [ ! $pathtoref ]; then
	printf "Missing arguments!\n"
	printf "$usage\n"
	exit 1
fi


#### To check each reference genomes for fasta and gff3 files, and make indexes if needed
### Not very great way (since it creates temporary unzipped fasta and gff3 files, even if not needed) but does the job quickly...
printf "\nMaking sure all reference directories are ready or if genome indexes need to be built\n"
ref_dir=$pathtoref/$ref
if [ -s $ref_dir/*.fa.gz ]; then
	fafile=$(ls $ref_dir/*.fa.gz)
	fafilename=${fafile##*/}
	printf "\nGzipped fasta file found in ${ref_dir}:\n $fafilename\n"
	pigz -p $threads -dc $fafile > $ref_dir/temp_${ref}.fa
	fasta=$ref_dir/temp_${ref}.fa
elif [ -s $ref_dir/*.fa ]; then
	fafile=$(ls $ref_dir/*.fa)
	fafilename=${fafile##*/}
	printf "\nUnzipped fasta file found in ${ref_dir}:\n $fafilename\n"
	fasta=$fafile
else
	printf "\nNo fasta file found in reference directory:\n ${ref_dir}\n"
	exit 1
fi

if [ -s $ref_dir/*.gff3.gz ]; then
	annofile=$(ls $ref_dir/*gff3.gz)
	annofilename=${annofile##*/}
	printf "\nGzipped annotation file found in ${ref_dir}:\n $annofilename\n"
	pigz -p $threads -dc $annofile > $ref_dir/temp_${ref}.gff3	
	gff=$ref_dir/temp_${ref}.gff3
elif [ -s $ref_dir/*.gff3 ]; then
	annofile=$(ls $ref_dir/*.gff3)
	annofilename=${annofile##*/}
	printf "\nUnzipped annotation file found in ${ref_dir}:\n $annofilename\n"
	gff=$annofile
else
	printf "\nNo annotation file (gff3) found in reference directory:\n ${ref_dir}\n"
	exit 1
fi

#### To create - if needed - the genome index, a bed file of all genes (for deeptools compatibility) and chrom.sizes file
if [ ! -s $ref_dir/chrom.sizes ]; then
	printf "\nMaking chrom.sizes file for $ref\n"
	samtools faidx $fasta
	cut -f1,2 ${fasta}.fai > $ref_dir/chrom.sizes
fi
	
if [ ! -s $datatype/tracks/${ref}_all_genes.bed ]; then
	printf "\nMaking a bed file with gene coordinates from $ref\n"
	awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,$9,".",$7}' $gff > $datatype/tracks/${ref}_all_genes.bed
fi

if [[ $datatype == "ChIP" ]]; then
	if [ ! -s $datatype/reports/summary_mapping_stats.txt ]; then
		printf "Line\tTissue\tMark\tRep\tTotal_reads\tPassing_filtering\tDeduplicated_reads\tProperly_mapped_reads\n" > $datatype/reports/summary_mapping_stats.txt
	fi
	if [ ! -e $ref_dir/*.bt2* ]; then
#### This return the following warning '[: too many arguments' when the index is already build (several bt2 files)
		printf "\nBuilding Bowtie2 index for $ref\n"
		bowtie2-build --threads $threads $fasta $ref_dir/$ref
	fi
elif [[ $datatype == "RNA" ]]; then
	if [ ! -s $datatype/reports/summary_mapping_stats.txt ]; then
		printf "Line\tTissue\tMark\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tOther\n" > $datatype/reports/summary_mapping_stats.txt
	fi
	gen_dir=$ref_dir/STAR_index
	if [ ! -d ${gen_dir} ]; then
		printf "Building STAR index directory for $ref\n"
		mkdir ${gen_dir}
		STAR --runThreadN $threads --runMode genomeGenerate --genomeDir ${gen_dir} --genomeFastaFiles $fasta --sjdbGTFfile $gff --sjdbGTFtagExonParentTranscript Parent
	fi
else
	printf "\nType of data unknown!\n"
	exit 1
fi

rm -f $ref_dir/temp*

printf "\nScript finished successfully!\n"
touch $datatype/chkpts/env_${ref}
