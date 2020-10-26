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
##### Main script for Maize code data analysis
##### 
##### sh MaizeCode.sh -f samplefile -t type -p path to genome reference
##### 	-f: samplefile (tab delimited text file with Line, Tissue, Sample, Rep, SequencingID, Path, PE or SE, Genome reference to map to)
##### 	-p: path to the folder containing all the different genome references (e.g. ~/data/Genomes/Zea_mays)
##### 	-h: help, returns usage
#####
##### Each reference genome should then be in a seperate folder (same name given in the last column of samplefile)
##### Each reference genome folder should contain a single fasta file and a single gff3 file (can be gzipped)
##### The gff3 files should have 'gene' in column 3 and exons should be linked by 'Parent' in column 9
##### The fasta and gff3 files should have the same chromosome names (i.e. 1 2 3... and 1 2 3... or Chr1 Chr2 Chr3... and Chr1 Chr2 Chr3...)
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
export mc_dir=$(dirname "$0")
echo "Running MaizeCode scripts from ${mc_dir} in working directory "${PWD}"

if [ $# -eq 0 ]; then
	printf "$usage\n"
	exit 1
fi

while getopts "f:t:p:h" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		p)	export pathtoref=${OPTARG};;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $samplefile ] || [ ! $pathtoref ]; then
	printf "Missing arguments!\n"
	printf "$usage\n"
	exit 1
fi

#############################################################################################
########################################### PART1 ###########################################
################################### Preparing enrivonment ###################################
#############################################################################################

#### To create the list of reference genomes needed and types of data to analyze
datatype_list=()
ref_list=()
datatype_ref_list=()
while read line tissue sample rep sampleID path paired ref
do
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
	esac
	datatype_list+=("$datatype")
	ref_list+=("$ref")
	data_ref_list+=("${datatype}_${ref}")
done < $samplefile
uniq_ref_list=($(printf "%s\n" "${ref_list[@]}" | sort -u))
uniq_datatype_list=($(printf "%s\n" "${datatype_list[@]}" | sort -u))


#### To check each reference genomes for fasta and gff3 files, and make indexes if needed
### Not very great way (since it creates temporary unzipped fasta and gff3 files, even if not needed) but does the job quickly...
printf "\nMaking sure all reference directories are ready or if genome indexes need to be built\n"
for ref in ${uniq_ref_list[@]}
do
	ref_dir=$pathtoref/$ref
	if [ -f $ref_dir/*.fa.gz ]; then
		fafile=$(ls $ref_dir/*.fa.gz)
		fafilename=${fafile##*/}
		printf "\nGzipped fasta file found in ${ref_dir}:\n $fafilename\n"
		pigz -p $threads -dc $fafile > $ref_dir/temp_${ref}.fa
		fasta=$ref_dir/temp_${ref}.fa
	elif [ -f $ref_dir/*.fa ]; then
		fafile=$(ls $ref_dir/*.fa)
		fafilename=${fafile##*/}
		printf "\nUnzipped fasta file found in ${ref_dir}:\n $fafilename\n"
		fasta=$fafile
	else
		printf "\nNo fasta file found in reference directory:\n ${ref_dir}\n"
		exit 1
	fi

	if [ -f $ref_dir/*.gff3.gz ]; then
		annofile=$(ls $ref_dir/*gff3.gz)
		annofilename=${annofile##*/}
		printf "\nGzipped annotation file found in ${ref_dir}:\n $annofilename\n"
		pigz -p $threads -dc $annofile > $ref_dir/temp_${ref}.gff3	
		gff=$ref_dir/temp_${ref}.gff3
	elif [ -f $ref_dir/*.gff3 ]; then
		annofile=$(ls $ref_dir/*.gff3)
		annofilename=${annofile##*/}
		printf "\nUnzipped annotation file found in ${ref_dir}:\n $annofilename\n"
		gff=$annofile
	else
		printf "\nNo annotation file (gff3) found in reference directory:\n ${ref_dir}\n"
		exit 1
	fi

	#### To create - if needed - the genome index, a bed file of all genes (for deeptools compatibility) and chrom.sizes file
	if [ ! -f $ref_dir/chrom.sizes ]; then
		printf "\nMaking chrom.sizes file for $ref\n"
		samtools faidx $fasta
		cut -f1,2 ${fasta}.fai > $ref_dir/chrom.sizes
	fi
	
	if [ ! -f $datatype/tracks/${ref}_all_genes.bed ]; then
		printf "\nMaking a bed file with gene coordinates from $ref\n"
		awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,$9,".",$7}' $gff > $datatype/tracks/${ref}_all_genes.bed
	fi

	for datatype in ${uniq_datatype_list[@]}
	do
		if [ ! -d ./$datatype ]; then
			mkdir ./$datatype
			mkdir ./$datatype/fastq
			mkdir ./$datatype/mapped
			mkdir ./$datatype/tracks
			mkdir ./$datatype/reports
		fi
		if [[ $datatype == "ChIP" ]] && [[ " ${data_ref_list[@]} " =~ " ${datatype}_${ref} " ]]; then
			if [ ! -f $datatype/reports/summary_mapping_stats.txt ]; then
				printf "Line\tTissue\tMark\tRep\tTotal_reads\tPassing_filtering\tDeduplicated_reads\tProperly_mapped_reads\n" > $datatype/reports/summary_mapping_stats.txt
			fi
			if [ ! -f $ref_dir/*.bt2* ]; then
#### This return the following warning '[: too many arguments' when the index is already build (several bt2 files)
				printf "\nBuilding Bowtie2 index for $ref\n"
				bowtie2-build --threads $threads $fasta $ref_dir/$ref
			fi
		elif [[ $datatype == "RNA" ]] && [[ " ${data_ref_list[@]} " =~ " ${datatype}_${ref} " ]]; then
			if [ ! -f $datatype/reports/summary_mapping_stats.txt ]; then
				printf "Line\tTissue\tMark\tRep\tTotal_reads\tPassing_filtering\tOther\n" > $datatype/reports/summary_mapping_stats.txt
			fi
			gen_dir=$ref_dir/STAR_index
			if [ ! -d ${gen_dir} ]; then
				printf "Building STAR index directory for $ref\n"
				mkdir ${gen_dir}
				STAR --runThreadN $threads --runMode genomeGenerate --genomeDir ${gen_dir} --genomeFastaFiles $fasta --sjdbGTFfile $gff --sjdbGTFtagExonParentTranscript Parent
			fi
		elif [[ ! " ${data_ref_list[@]} " =~ " ${datatype}_${ref} " ]]; then
			### Combination datetype * ref does not exist in the sample file, moving on
			:
		else
			printf "\nType of data unknown!\n"
			exit 1
		fi
	done
done
rm -f $ref_dir/temp*


#############################################################################################
########################################### PART2 ###########################################
##################################### Getting samples #######################################
#############################################################################################


while read line tissue sample rep sampleID path paired ref
do
	ref_dir=$pathtoref/$ref
	name=${line}_${tissue}_${sample}_${rep}
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
	esac

	if [ ! -f ./$datatype/fastq/${name}*.fastq.gz ]; then
#### This return the following warning '[: too many arguments' when several files are there (PE data)
		if [[ $paired == "PE" ]]; then
			printf "\nCopying PE fastq for $name ($sampleID in $path)\n"
			cp $path/${sampleID}*R1*fastq.gz ./$datatype/fastq/${name}_R1.fastq.gz
			cp $path/${sampleID}*R2*fastq.gz ./$datatype/fastq/${name}_R2.fastq.gz
		elif [[ $paired == "SE" ]]; then
			printf "\nCopying SE fastq for $name ($sampleID in $path)\n"
			cp $path/${sampleID}*fastq.gz ./$datatype/fastq/${name}.fastq.gz
		else
			printf "\nData format unknown: paired-end (PE) or single-end (SE)?\n"
			exit 1
		fi
	fi
	printf "\nRunning $datatype worker script for $name\n"
	cd $datatype
	qsub -N ${name} -o ${name}.log ~/data/Scripts/MaizeCode_${datatype}_sample.sh -d $ref_dir -l $line -t $tissue -m $sample -r $rep -p $paired
	cd ..
done < $samplefile


############################################################################################
########################################## MISC ############################################
############################################################################################

# # ##### To make the samplefile (e.g. for B73 endosperm)

# # printf "B73\tendosperm\tInput\tRep1\tS01\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me1\tRep1\tS02\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me3\tRep1\tS03\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K27ac\tRep1\tS04\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tInput\tRep2\tS05\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me1\tRep2\tS06\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me3\tRep2\tS07\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K27ac\tRep2\tS08\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\n" > B73_endosperm_samplefile.txt

############################################################################################

printf "\nScript finished successfully!\n"
