#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=2G
#$ -l tmp_free=20G
#$ -o env.log
#$ -j y
#$ -N env

usage="
##### Script for checking environment and building indexes for Maize code
##### 
##### sh MaizeCode_check_environment.sh -p path to genome reference -r ref -d datatype
##### 	-p: path to the folder containing all the different genome references (e.g. ~/data/Genomes/Zea_mays)
##### 	-r: genome reference to use (e.g. B73_v4) 
#####	-d: type of data [ ChIP | RNA | TF | shRNA | mC ]
##### 	-h: help, returns usage
#####
##### The reference genome folder should contain a single fasta file (.fa or .fasta), a single GFF file (.gff [or .gff*]) and a single GTF file (can be gzipped)
##### The GFF file should have 'gene' in column 3
##### The GTF file can be produced using cufflinks with gffread -T <gff_file> -o <gtf_file>
##### The fasta, GFF and GTF files should have the same chromosome names (i.e. 1 2 3... and 1 2 3... or Chr1 Chr2 Chr3... and Chr1 Chr2 Chr3...)
##### The type of data (ChIP or RNA) will define whether to build bowtie2 or STAR indexes
#####
##### This script check if the indexes and other files need to be created and creates them if they do
#####
##### Requirements: pigz, samtools, Bowtie2 for ChIP data, STAR for RNA data, bismark for mC data
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

while getopts "p:r:d:h" opt; do
	case ${opt} in
		h) 	printf "${usage}\n"
			exit 0;;
		p)	export pathtoref=${OPTARG};;
		r)	export ref=${OPTARG};;
		d)	export datatype=${OPTARG};;
		*)	printf "${usage}\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! ${ref} ] || [ ! ${pathtoref} ] || [ ! ${datatype} ]; then
	printf "Missing arguments!\n"
	printf "${usage}\n"
	exit 1
fi


#### To check each reference genomes for fasta, GFF3 and GTF files, and make indexes if needed
printf "\nMaking sure all reference directories are ready or if genome indexes need to be built\n"
ref_dir=${pathtoref}/${ref}
if [ -s ${ref_dir}/*.fa.gz ]; then
	fa_file=$(ls ${ref_dir}/*.fa.gz)
	fa_filename=${fa_file##*/}
	printf "\nGzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
	pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${datatype}_${ref}.fa
	fasta=${ref_dir}/temp_${datatype}_${ref}.fa
elif [ -s ${ref_dir}/*.fa ]; then
	fa_file=$(ls ${ref_dir}/*.fa)
	fa_filename=${fa_file##*/}
	printf "\nUnzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
	fasta=${fa_file}
elif [ -s ${ref_dir}/*.fasta.gz ]; then
	fa_file=$(ls ${ref_dir}/*.fasta.gz)
	fa_filename=${fa_file##*/}
	printf "\nGzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
	pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${datatype}_${ref}.fa
	fasta=${ref_dir}/temp_${datatype}_${ref}.fa
elif [ -s ${ref_dir}/*.fasta ]; then
	fa_file=$(ls ${ref_dir}/*.fasta)
	fa_filename=${fa_file##*/}
	printf "\nUnzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
	fasta=${fa_file}
else
	printf "\nNo fasta file found in reference directory:\n ${ref_dir}\n"
	exit 1
fi

if [ -s ${ref_dir}/*.gff*.gz ]; then
	gff_file=$(ls ${ref_dir}/*gff*.gz)
	gff_filename=${gff_file##*/}
	printf "\nGzipped GFF annotation file found in ${ref_dir}:\n ${gff_filename}\n"
	pigz -p ${threads} -dc ${gff_file} > ${ref_dir}/temp_${datatype}_${ref}.gff	
	gff=${ref_dir}/temp_${datatype}_${ref}.gff
elif [ -s ${ref_dir}/*.gff* ]; then
	gff_file=$(ls ${ref_dir}/*.gff*)
	gff_filename=${gff_file##*/}
	printf "\nUnzipped GFF annotation file found in ${ref_dir}:\n ${gff_filename}\n"
	gff=${gff_file}
else
	printf "\nNo gff annotation file found in reference directory:\n ${ref_dir}\n"
	exit 1
fi

if [ -s ${ref_dir}/*.gtf.gz ]; then
	gtf_file=$(ls ${ref_dir}/*gtf.gz)
	gtf_filename=${gtf_file##*/}
	printf "\nGzipped GTF annotation file found in ${ref_dir}:\n ${gtf_filename}\n"
	pigz -p ${threads} -dc ${gtf_file} > ${ref_dir}/temp_${datatype}_${ref}.gtf	
	gtf=${ref_dir}/temp_${datatype}_${ref}.gtf
elif [ -s ${ref_dir}/*.gtf ]; then
	gtf_file=$(ls ${ref_dir}/*.gtf)
	gtf_filename=${gtf_file##*/}
	printf "\nUnzipped GTF annotation file found in ${ref_dir}:\n ${gtf_filename}\n"
	gtf=${gtf_file}
else
	printf "\nNo GTF annotation file found in reference directory:\n ${ref_dir}\n"
	exit 1
fi

#### To create - if needed - the genome index, a bed file of all genes (for future analysis) and a chrom.sizes file
if [ ! -s ${ref_dir}/chrom.sizes ]; then
	printf "\nMaking chrom.sizes file for ${ref}\n"
	samtools faidx ${fasta}
	cut -f1,2 ${fasta}.fai > ${ref_dir}/chrom.sizes
fi

if [ ! -s ${datatype}/tracks/${ref}_all_genes.bed ]; then
	printf "\nMaking a bed file with gene coordinates from ${ref}\n"
	awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,$9,".",$7}' ${gff} | bedtools sort -g ${ref_dir}/chrom.sizes > ${datatype}/tracks/${ref}_protein_coding_genes.bed
	awk -v OFS="\t" '$3~"gene" {print $1,$4-1,$5,$9,".",$7}' ${gff} | bedtools sort -g ${ref_dir}/chrom.sizes > ${datatype}/tracks/${ref}_all_genes.bed
fi

if [[ ${datatype} == "ChIP" ]]; then
	if [ ! -s ${datatype}/reports/summary_mapping_stats.txt ]; then
		printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > ${datatype}/reports/summary_mapping_stats.txt
	fi
	if ls ${ref_dir}/*.bt2* 1> /dev/null 2>&1; then
		printf "\nBowtie2 index already exists for ${ref} in ${ref_dir}\n"
	elif grep -q "Building Bowtie2 index" TF/logs/env_${ref}.log; then
		printf "\nBowtie2 index already being prepared for TF for ${ref} in ${ref_dir}\n"
	else
		printf "\nBuilding Bowtie2 index for $ref\n"
		bowtie2-build --threads ${threads} ${fasta} ${ref_dir}/${ref}
	fi
elif [[ ${datatype} == "RNA" ]]; then
	if [ ! -s ${datatype}/reports/summary_mapping_stats.txt ]; then
		printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > ${datatype}/reports/summary_mapping_stats.txt
	fi
	gen_dir=${ref_dir}/STAR_index
	if [ ! -d ${gen_dir} ]; then
		printf "\nBuilding STAR index directory for ${ref}\n"
		mkdir ${gen_dir}
		STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir ${gen_dir} --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf}
	else
		printf "\nSTAR index already exists for ${ref} in ${ref_dir}\n"
	fi
elif [[ ${datatype} == "TF" ]]; then
	if [ ! -s ${datatype}/reports/summary_mapping_stats.txt ]; then
		printf "Line\tSample\tMark\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > ${datatype}/reports/summary_mapping_stats.txt
	fi
	if ls ${ref_dir}/*.bt2* 1> /dev/null 2>&1; then
		printf "\nBowtie2 index already exists for ${ref} in ${ref_dir}\n"
	elif grep -q "Building Bowtie2 index" ChIP/logs/env_${ref}.log; then
		printf "\nBowtie2 index already being prepared for ChIP for ${ref} in ${ref_dir}\n"
	else
		printf "\nBuilding Bowtie2 index for ${ref}\n"
		bowtie2-build --threads ${threads} ${fasta} ${ref_dir}/${ref}
	fi
	#### This step will need to be automatized to potentially change which line/organism to have the annotations to mask.
	#### Another option would be to have people create it independantly with Repeat Masker, potentially giving a script or just help with documentation.
	if [[ ${ref} == "B73_v4" ]] && [ ! -s ${datatype}/tracks/B73_v4_masked_regions.bed ]; then
		printf "\nCopying bed files of masked regions to filter peaks\n"
		cp /grid/martienssen/data_norepl/dropbox/maizecode/Repeat_Masker/B73_v4_masked_regions.bed ${datatype}/tracks/
	fi
elif [[ ${datatype} == "mC" ]]; then
	if [ ! -d ${datatype}/methylcall ]; then
		mkdir ${datatype}/methylcall
	fi
	if [ ! -s ${datatype}/reports/summary_mapping_stats.txt ]; then
		printf "Line\tTissue\tRep\tTotal_reads\tAll_mapped_reads\tUniquely_mapped_reads\tPercentage_covered\tPercentage_covered_min3reads\tAverage_coverage_all\tAverage_coverage_covered\tNon_conversion_rate(Pt/Lambda)\n" > ${datatype}/reports/summary_mapping_stats.txt
	fi
	if [ ! -d ${ref_dir}/Bisulfite_Genome ]; then
		printf "\nBuilding index...\n"
		bismark_genome_preparation --bowtie2 --genomic_composition ${ref_dir}
		printf "\nBuilding index finished!\n"
	fi
elif [[ ${datatype} == "shRNA" ]]; then
	if ls ${ref_dir}/*.bt2* 1> /dev/null 2>&1; then
		printf "\nBowtie2 index already exists for ${ref} in ${ref_dir}\n"
	elif grep -q "Building Bowtie2 index" ChIP/logs/env_${ref}.log; then
		printf "\nBowtie2 index already being prepared for ChIP for ${ref} in ${ref_dir}\n"
	elif grep -q "Building Bowtie2 index" TF/logs/env_${ref}.log; then
		printf "\nBowtie2 index already being prepared for TF for ${ref} in ${ref_dir}\n"
	else
		printf "\nBuilding Bowtie2 index for ${ref}\n"
		bowtie2-build --threads ${threads} ${fasta} ${ref_dir}/${ref}
	fi
	#### This step will need to be automatized to potentially change which organism to map to and build it from scratch. 
	#### Another option would be to have people create it independantly, potentially giving a script or just the help documentation (see Help_RFam.docx in the data folder).
	printf "\nCopying structural RNA index files to filter reads\n"
	cp -r /grid/martienssen/data_nlsas_norepl/dropbox/maizecode/structural_RNA ${datatype}/
else
	printf "\nType of data unknown!\n"
	exit 1
fi

rm -f ${ref_dir}/temp_${datatype}*

printf "\nScript finished successfully!\n"
touch ${datatype}/chkpts/env_${ref}

