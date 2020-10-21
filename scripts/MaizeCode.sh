#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 8
#$ -l m_mem_free=10G
#$ -l tmp_free=50G
#$ -o maizecode.log
#$ -j y
#$ -N maizecode

usage="
##### Master script for Maize code data analysis
##### 
##### sh MaizeCode.sh -f samplefile -t type -p path to genome reference -r name of reference genome
##### 	-f: samplefile (tab delimited text file with Line, Tissue, Sample, Rep, Path, PE or SE)
##### 	-t: data type (ChIP, RNA or RAMPAGE)
##### 	-p: path to the folder containing all the different genome references (e.g. ~/data/Genomes/Zea_mays)
##### 	-r: reference genome (line ID + ref; e.g. B73_v4)
##### 	-h: help, returns usage
#####
##### Each reference genome should then be in a seperate folder (named like argument #4) with .fa and .gff3 files (can be gzipped)
#####
##### This script creates the folders needed, copies the fastq files, 
##### prepares the genome index and chrom.sizes file if they are not done for this type of data
##### and send each sample into MaizeCode_type.in.arg2_sample.sh
#####
##### Requirements: pigz, samtools, Bowtie2 for ChIP data, STAR for RNA data, xxx for RAMPAGE data
##### ! Need to check or change the location of your Scripts !
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

while getopts "f:t:p:r:h" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		t)	export type=${OPTARG};;
		p)	export pathtoref=${OPTARG};;
		r)	export ref=${OPTARG};;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $samplefile ] || [ ! $type ] || [ ! $pathtoref ] || [ ! $ref ]; then
	printf "Missing arguments!\n"
	printf "$usage\n"
	exit 1
fi

if [ ! -d ./$type ]; then
	mkdir ./$type
	mkdir ./$type/fastq
	mkdir ./$type/mapped
	mkdir ./$type/deeptools
	mkdir ./$type/reports
fi

ref_dir=$pathtoref/$ref
printf "\nReference directory containing genome files: $ref_dir\n"

#############################################################################################
########################################### PART1 ###########################################
################################### Preparing enrivonment ###################################
#############################################################################################

### Not very great way (since it creates temporary unzipped fasta and gff3 files, even if not needed) but does the job quickly...

printf "\nMaking sure reference directory is ready or if genome indexes need to be built\n"

if [ -f $ref_dir/*.fa.gz ]; then
	fafile=$(ls $ref_dir/*.fa.gz)
	printf "\nGzipped fasta file found: $fafile\n"
	pigz -p $threads -dc $fafile > $ref_dir/temp_${ref}.fa
	fasta=$ref_dir/temp_${ref}.fa
elif [ -f $ref_dir/*.fa ]; then
	fafile=$(ls $ref_dir/*.fa)
	printf "\nUnzipped fasta file found: $fafile\n"
	fasta=$fafile
else
	printf "\nNo fasta file found in reference directory\n"
	exit 1
fi

if [ -f $ref_dir/*.gff3.gz ]; then
	annofile=$(ls $ref_dir/*gff3.gz)
	printf "\nGzipped annotation file found: $annofile\n"
	pigz -p $threads -dc $annofile > $ref_dir/temp_${ref}.gff3	
	gff=$ref_dir/temp_${ref}.gff3
elif [ -f $ref_dir/*.gff3 ]; then
	annofile=$(ls $ref_dir/*.gff3)
	printf "\nUnzipped annotation file found: $annofile\n"
	gff=$annofile
else
	printf "\nNo annotation file (gff3) found in reference directory\n"
	exit 1
fi

#### To create - if needed - the genome index, a bed file of all genes (for deeptools compatibility) and chrom.sizes file
if [ ! -f $ref_dir/chrom.sizes ]; then
	printf "\nMaking chrom.sizes file for $ref\n"
	samtools faidx $fasta
	cut -f1,2 ${fasta}.fai > $ref_dir/chrom.sizes
fi
	
if [ ! -f $type/deeptools/all_genes.bed ]; then
	awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,$9,".",$7}' $gff > $type/deeptools/${ref}_all_genes.bed
fi

if [[ $type == "ChIP" ]]; then
	if [ ! -f $ref_dir/*.bt2* ]; then
#### This return the following warning 'line 111: [: too many arguments when the index is already build'
		printf "\nBuilding Bowtie2 index for $ref\n"
		bowtie2-build --threads $threads $fasta $ref_dir/$ref
	fi
elif [[ $type == "RNA" ]]; then
	gen_dir=$ref_dir/STAR_index
	if [ ! -d ${gen_dir} ]; then
		printf "Building STAR index directory for $ref\n"
		mkdir ${gen_dir}
		STAR --runThreadN $threads --runMode genomeGenerate --genomeDir ${gen_dir} --genomeFastaFiles $fasta --sjdbGTFfile $gff --sjdbGTFtagExonParentTranscript Parent
	fi
elif [[ $type == "RAMPAGE" ]]; then
	printf "\nRAMPAGE mode not yet activated.. Sorry!\n"
else
	printf "\nType of data unknown!\n"
	exit 1
fi

rm -f $ref_dir/temp*


#############################################################################################
########################################### PART2 ###########################################
##################################### Getting samples #######################################
#############################################################################################


if [ ! -f $type/reports/summary_mapping_stats.txt ]; then
	printf "Line\tTissue\tMark\tRep\tTotal_reads\tPassing_filtering\tDeduplicated_reads\tProperly_mapped_reads\n" > $type/reports/summary_mapping_stats.txt
fi

while read line tissue sample rep sampleID path paired
do
	name=${line}_${tissue}_${sample}_${rep}
	if [ ! -f ./$type/fastq/${name}*.fastq.gz ]; then
		if [[ $paired == "PE" ]]; then
			printf "\nCopying PE fastq for $name ($sampleID in $path)\n"
			cp $path/${sampleID}*R1*fastq.gz ./$type/fastq/${name}_R1.fastq.gz
			cp $path/${sampleID}*R2*fastq.gz ./$type/fastq/${name}_R2.fastq.gz
		elif [[ $paired == "SE" ]]; then
			printf "\nCopying SE fastq for $name ($sampleID in $path)\n"
			cp $path/${sampleID}*fastq.gz ./$type/fastq/${name}.fastq.gz
		else
			printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
			exit 1
		fi
	fi
	printf "\nRunning $type worker script for $name\n"
	cd $type
	qsub -N ${name} -o ${name}.log ~/data/Scripts/MaizeCode_${type}_sample.sh -d $ref_dir -l $line -t $tissue -m $sample -r $rep -p $paired
	cd ..
done < $samplefile


############################################################################################
########################################## MISC ############################################
############################################################################################

# # ##### To make the samplefile (e.g. for B73 endosperm)

# # printf "B73\tendosperm\tInput\tRep1\tS01\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\nB73\tendosperm\tH3K4me1\tRep1\tS02\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\nB73\tendosperm\tH3K4me3\tRep1\tS03\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\nB73\tendosperm\tH3K27ac\tRep1\tS04\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\nB73\tendosperm\tInput\tRep2\tS05\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\nB73\tendosperm\tH3K4me1\tRep2\tS06\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\nB73\tendosperm\tH3K4me3\tRep2\tS07\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\nB73\tendosperm\tH3K27ac\tRep2\tS08\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\n" > B73_endosperm_samplefile.txt

############################################################################################

printf "\nScript finished successfully!\n"
