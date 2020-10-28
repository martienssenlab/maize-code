#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=2G
#$ -l tmp_free=10G
#$ -o maizecode.log
#$ -j y
#$ -N maizecode

usage="
##### Main script for Maize code data analysis
##### 
##### sh MaizeCode.sh -f samplefile -p path to genome reference [-s]
##### 	-f: samplefile
##### 	-p: path to the folder containing all the different genome references (e.g. ~/data/Genomes/Zea_mays)
#####	-s: if set, the analysis does not proceed (default=not set, keep going with the analysis over all the samples in the samplefile)
##### 	-h: help, returns usage
#####
##### The samplefile should be a tab delimited text file with 8 columns:
##### col #1: Line (e.g. B73)
##### col #2: Tissue (e.g endosperm) 
##### col #3: Sample (e.g. H3K4me3 or 'Input' for ChIP, (sh)RNA or RAMPAGE for RNA). Histone marks must start with capital 'H'. 'Input', RNA and RAMPAGE are case-sensitive
##### col #4: Replicate ID (e.g. Rep1 or Rep2). If only one replicate in the experiment, put 'Rep1'.
##### col #5: SequencingID (e.g. S01). Unique identifier for the sample in the sequencing folder
##### col #6: Path to the fastq files (e.g. /seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846)
##### col #7: If data is paired-end or single-end [PE | SE]. 
##### col #8: Name of the genome reference to map. Each genome reference should have a unique folder that contains a single fasta file and a single gff3 file (can be gzipped).
##### The gff3 files should have 'gene' in column 3 and exons should be linked by 'Parent' in column 9
##### The fasta and gff3 files should have the same chromosome names (i.e. 1 2 3... and 1 2 3... or Chr1 Chr2 Chr3... and Chr1 Chr2 Chr3...)
##### For cleaner naming purposes, use '_samplefile.txt' as suffix
#####
##### This script creates the folders needed,
##### prepares the genome index and chrom.sizes file if they are not done for this type of data
##### and send each sample into a datatype-specific script (MaizeCode_$datatype_sample.sh) if they have not been mapped before
#####
##### Requirements: pigz, samtools, Bowtie2 for ChIP data, STAR for RNA data, xxx for RAMPAGE data
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS
# # export mc_dir=$(dirname "$0")
export mc_dir="${HOME}/data/Scripts/MaizeCode"
printf "\nRunning MaizeCode scripts from ${mc_dir} in working directory ${PWD}\n"

if [ $# -eq 0 ]; then
	printf "$usage\n"
	exit 1
fi

while getopts "f:p:sh" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		p)	export pathtoref=${OPTARG};;
		s)	printf "\nOption not to perform analysis selected\n"
			export keepgoing="STOP";;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $samplefile ]; then
	printf "Samplefile missing!\n"
	printf "$usage\n"
	exit 1
fi
if [ ! $pathtoref ]; then
	printf "Path to reference genome folders missing!\n"
	printf "$usage\n"
	exit 1
fi

#############################################################################################
########################################### PART1 ###########################################
################################### Preparing enrivonments ##################################
#############################################################################################

if [ ! -d ./chkpts ]; then
	mkdir ./chkpts
fi

#### To prepare the lists of samples, reference genomes and type of data

newdatatype_list=()
newref_list=()
datatype_ref_list=()
new_env=0
new_sample=0
while read line tissue sample rep sampleID path paired ref
do
	name=${line}_${tissue}_${sample}_${rep}
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
		*) datatype="unknown";;
	esac
	if [[ "$datatype" == "unknown" ]]; then
		printf "Type of data unknown!\n"
		printf "$usage\n"
		exit 1
	fi
	if [ ! -e $datatype/chkpts/${name}_${ref} ]; then
		new_sample=1
	fi
	if [ ! -e $datatype/chkpts/env_${ref} ]; then
		new_env=1
		newdatatype_list+=("$datatype")
		newref_list+=("$ref")
		data_ref_list+=("${datatype}_${ref}")
	fi
done < $samplefile


#### To check if new environments need to be prepapred
if [[ ${new_env} == 0 ]]; then
	printf "\nAll environments are ready for mapping\n"
else
	uniq_newref_list=($(printf "%s\n" "${newref_list[@]}" | sort -u))
	uniq_newdatatype_list=($(printf "%s\n" "${newdatatype_list[@]}" | sort -u))

	check_list=()
	pids=()
	for ref in ${uniq_newref_list[@]}
	do
		for datatype in ${uniq_newdatatype_list[@]}
		do
			if [[ " ${data_ref_list[@]} " =~ " ${datatype}_${ref} " ]]; then
				check_list+=("$datatype/chkpts/env_${ref}")
				if [ ! -d ./$datatype ]; then
					mkdir ./$datatype
					mkdir ./$datatype/fastq
					mkdir ./$datatype/mapped
					mkdir ./$datatype/tracks
					mkdir ./$datatype/reports
					mkdir ./$datatype/logs
					mkdir ./$datatype/chkpts
				fi
				printf "\nPreparing environment of ${ref} genome for ${datatype} data\n"
				qsub -sync y -N env_${ref}_${datatype} -o $datatype/logs/env_${ref}.log ${mc_dir}/MaizeCode_check_environment.sh -p $pathtoref -r $ref -d $datatype &
				pids+=("$!")
			elif [[ ! " ${data_ref_list[@]} " =~ " ${datatype}_${ref} " ]]; then
			### Combination datetype * ref does not exist in the sample file, moving on
				:
			elif [ ! -d $pathtoref/$ref ]; then
				printf "\nNo $ref folder found in $pathtoref\n"
				printf "$usage\n"
				exit 1
			else
				printf "\nError!\n"
				printf "$usage\n"
				exit 1
			fi	
		done
	done
	#### Wait for the check_environment scripts to finish
	printf "\nWaiting for the environments to be prepared...\n"
	wait $pids
	#### Check if the environment are good or if an error occurred
	for check in ${check_list[@]}
	do
		datatype=${check%%/*}
		ref=${check##*/}
		if [ ! -e $check ]; then
			printf "\nProblem while making environment for $datatype with $ref genome!\nCheck log: $datatype/logs/env_${ref}.log\n"
			exit 1
		else
			printf "\nEnvironment for $datatype with $ref genome now good to go\n"
		fi
	done
fi


#############################################################################################
########################################### PART2 ###########################################
############################# Getting samples, mapping and QC ###############################
#############################################################################################

tmp1=${samplefile%%_samplefile*}
analysisfile="${tmp1}_analysis_samplefile.txt"

if [ -s ${analysisfile} ]; then
	rm -f ${analysisfile}
fi

check_list=()
checkname_list=()
checkdatatype_list=()
ref_list=()
pids=()
while read line tissue sample rep sampleID path paired ref
do
	ref_dir=$pathtoref/$ref
	name=${line}_${tissue}_${sample}_${rep}
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
	esac
	check=$datatype/chkpts/${name}_${ref}
	checkname_list+=("$name")
	checkdatatype_list+=("$datatype")
	check_list+=("$check")
	ref_list+=("$ref")
	if [ -e $check ]; then
		printf "Sample $name has already been mapped to $ref genome\n"
	else
		if [ ! -e ./$datatype/fastq/${name}*.fastq.gz ]; then
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
		printf "\nRunning $datatype mapping script for $name on $ref genome\n"
		cd $datatype
		qsub -sync y -N ${name} -o logs/${name}.log ${mc_dir}/MaizeCode_${datatype}_sample.sh -d $ref_dir -l $line -t $tissue -m $sample -r $rep -p $paired &
		pids+=("$!")
		cd ..
	fi
	if [[ "${sample}" != "Input" ]]; then
		printf "${line}\t${tissue}\t${sample}\t${rep}\t${paired}\n" >> $analysisfile
	fi
done < $samplefile

if [[ ${new_sample} != 0 ]]; then
#### Wait for the mapping sample scripts to finish
	printf "\nWaiting for the samples to be mapped...\n"
	wait $pids

	i=0
	for check in ${check_list[@]}
	do
		if [ ! -e $check ]; then
			printf "\nProblem during mapping of ${checkname_list[i]}!\nCheck log: ${checkdatatype_list[i]}/logs/${checkname_list[i]}.log"
			exit 1
		else
			printf "\nSample ${checkname_list[i]} successfully mapped.\n"
		fi
		i=$((i+1))
	done
fi

#############################################################################################
########################################### PART3 ###########################################
############################# Starting analysis if not stopped ##############################
#############################################################################################

if [[ "$keepgoing" == "STOP" ]]; then
	printf "\nScript finished successfully without analysis\n"		
	exit 0
fi	

tmp1=${samplefile##*/}
samplename=${tmp1%%_samplefile*}

uniq_ref_list=($(printf "%s\n" "${ref_list[@]}" | sort -u))

pids=()
if [ ${#uniq_ref_list[@]} -eq 1 ]; then
	regionfile="${checkdatatype_list[0]}/tracks/${uniq_ref_list[0]}_all_genes.bed"
	tmp2=${regionfile##*/}
	regionname=${tmp2%.*}
	analysisname="${samplename}_on_${regionname}"
	check="combined/chkpts/${analysisname}"
	printf "\nPerforming the complete analysis using $regionname as region file\n"
	qsub -sync y -N maizecodeanalysis -o maizecode.log ${mc_dir}/MaizeCode_analysis.sh -f $analysisfile -r $regionfile &
	pids+=("$!")
else	
	analysisname="${samplename}_no_region"
	check="chkpts/${analysisname}"
	printf "\nToo many references for a combined analysis to be performed\nYou will need to run the MaizeCode_analysis.sh script once that this one is finished\nUsing a specific regionfile for combined analysis!\n"
	qsub -sync y -N maizecodeanalysis -o maizecode.log ${mc_dir}/MaizeCode_analysis.sh -f $analysisfile &
	pids+=("$!")
fi

wait $pids
if [ ! -e ${check} ]; then
	printf "\nProblem during the analysis. Check maizecode.log\n"
	exit 1
else
	printf "\nMaizeCode script finished successfully!\nCheck out the stats and plots!\n"
fi

############################################################################################
########################################## MISC ############################################
############################################################################################

# # ##### To make the samplefile (e.g. for B73 endosperm)

# # printf "B73\tendosperm\tInput\tRep1\tS01\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me1\tRep1\tS02\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me3\tRep1\tS03\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K27ac\tRep1\tS04\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tInput\tRep2\tS05\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me1\tRep2\tS06\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me3\tRep2\tS07\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K27ac\tRep2\tS08\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\n" > B73_endosperm_samplefile.txt

############################################################################################

