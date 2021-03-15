#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=2G
#$ -o maizecode.log
#$ -j y
#$ -N maizecode

usage="
##### Main script for Maize code data analysis
##### 
##### sh MaizeCode.sh -f <samplefile> -p <path to genome reference> [-s] [-c] [-h]
##### 	-f: samplefile
##### 	-p: path to the folder containing all the different genome references (e.g. ~/data/Genomes/Zea_mays)
#####	-s: if set, the whole analysis does not proceed (default=not set, keep going with the analysis over all the samples in the samplefile)
#####	-c: if set, only single samples analysis proceeds, not grouped analysis per line (default=not set, keep going with the complete analysis)
##### 	-h: help, returns usage
#####
##### The samplefile should be a tab-delimited text file with 8 columns:
##### col #1: Line (e.g. B73)
##### col #2: Tissue (e.g endosperm) 
##### col #3: Sample (e.g. H3K4me3 or 'Input' for ChIP, (sh)RNA or RAMPAGE for RNA). Histone marks must start with capital 'H'. 'Input', 'RNA' and 'RAMPAGE' are case-sensitive.
##### col #4: Replicate ID [Rep1 | Rep2]
##### col #5: SequencingID (e.g. S01). Unique identifier for the name of the sample in the raw sequencing folder which path is given in the next column. If downloading from SRA, put the SRR ID here.
##### col #6: Path to the fastq files (e.g. /seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846). If downloading from SRA, put 'SRA'.
##### col #7: If data is paired-end or single-end [PE | SE]. 
##### col #8: Name of the genome reference to map (e.g. B73_v4). Each genome reference should have a unique folder that contains a single fasta file and a single gff3 file (can be gzipped).
##### The gff3 files should have 'gene' in column 3 and exons should be linked by 'Parent' in column 9
##### The fasta and gff3 files should have the same chromosome names (i.e. 1 2 3... and 1 2 3... or Chr1 Chr2 Chr3... and Chr1 Chr2 Chr3...)
##### For cleaner naming purposes, use '_samplefile.txt' as suffix
#####
##### This script creates the folders needed,
##### prepares the genome index and chrom.sizes file if they are not done for this type of data,
##### sends each sample into a datatype-specific script for mapping (MaizeCode_$datatype_sample.sh) if they have not been mapped before
##### starts the script (MaizeCode_analysis.sh) for the analysis on all the samples in the samplefile (if -s is not set), 
##### It uses all the genes of the reference genome (if all samples are mapping to the same reference) as a region file or a combined analysis,
##### or only proceed with single sample analysis if different references are used. In the latter case, MaizeCode_analysis.sh script will need to be run independantly with the regionfile of your choice.
#####
##### Requirements for the mapping pipeline: pigz, samtools, fastQC, Cutadapt, Bowtie2 for ChIP data, STAR for RNA data, R and R packages: ggplot2,dplyr,tidyr,RColorBrewer,cowplot), parallel-fastq-dump (if downloading from SRA)
##### Additional requirements for the analysis pipeline: bedtools, deeptools, macs2, idr, R packages: UpSetR for ChIP and RNA data, and limma,edgeR,stringr,gplots for RNAseq data
"

set -e -o pipefail

printf "\n\n"
date
startdate=`date +%s`
printf "\n"

export threads=$NSLOTS
export mc_dir="${PWD}/scripts/"
printf "\nRunning MaizeCode.sh script in working directory ${PWD}\n"

if [ $# -eq 0 ]; then
	printf "$usage\n"
	exit 1
fi

while getopts "f:p:sch" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		p)	export pathtoref=${OPTARG};;
		s)	printf "\nOption not to perform analysis selected\n"
			export keepgoing="STOP";;
		c)	printf "\nOption not to perform combined analysis selected\n"
			export wholeanalysis="STOP";;
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
datat_ref_list=()
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
			### Combination datatype * ref does not exist in the sample file, moving on
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
	wait ${pids[*]}
	#### Check if the environment are good or if an error occurred
	for check in ${check_list[@]}
	do
		datatype=${check%%/*}
		ref=${check##*/env_}
		if [ ! -e $check ]; then
			printf "\nProblem while making environment for $datatype with $ref genome!\nCheck log: $datatype/logs/env_${ref}.log\n"
			exit 1
		else
			printf "\nEnvironment for $datatype with $ref genome now good to go\n"
			if [ -e ${pathtoref}/${ref}/temp_${ref}.fa ]; then
				rm -f ${pathtoref}/${ref}/temp_${ref}.fa
			fi
		fi
	done
fi

if [ ! -d ./combined ]; then
	mkdir ./combined
	mkdir ./combined/peaks
	mkdir ./combined/DEG
	mkdir ./combined/TSS
	mkdir ./combined/reports
	mkdir ./combined/matrix
	mkdir ./combined/plots
	mkdir ./combined/chkpts
	mkdir ./combined/logs
fi

#############################################################################################
########################################### PART2 ###########################################
############################# Getting samples, mapping and QC ###############################
#############################################################################################

tmp1=${samplefile%%_samplefile*}
analysisfile="${tmp1}_analysis_samplefile.txt"
samplename=${tmp1##*/}

if [ -s ${analysisfile} ]; then
	rm -f ${analysisfile}
fi

check_list=()
checkname_list=()
checkdatatype_list=()
ref_list=()
sample_list=()
pids=()
while read line tissue sample rep sampleID path paired ref
do
	ref_dir=$pathtoref/$ref
	shortname=${line}_${tissue}_${sample}
	name=${line}_${tissue}_${sample}_${rep}
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
	esac
	check=$datatype/chkpts/${name}_${ref}
	ref_list+=("$ref")
	if [ -e $check ]; then
		printf "Sample $name has already been mapped to $ref genome\n"
	else
		checkname_list+=("$name")
		checkdatatype_list+=("$datatype")
		check_list+=("$check")
		if ls ./$datatype/fastq/trimmed_${name}*.fastq.gz 1> /dev/null 2>&1; then
			printf "\nTrimmed fastq file(s) for ${name} already exist\n"
			export step="done"
		elif ls ./$datatype/fastq/${name}*.fastq.gz 1> /dev/null 2>&1; then
			printf "\nFastq file(s) for ${name} already exist\n"
			export step="trim"
		else
			printf "\nNew sample ${name} to be copied/downloaded\n"
			export step="download"
		fi
		printf "\nRunning $datatype mapping script for $name on $ref genome\n"
		cd $datatype
		qsub -sync y -N ${name} -o logs/${name}.log ${mc_dir}/MaizeCode_${datatype}_sample.sh -d $ref_dir -l $line -t $tissue -m $sample -r $rep -i $sampleID -f $path -p $paired -s $step &
		pids+=("$!")
		cd ..
	fi
	if [[ ! "${sample_list[@]}" =~ "${shortname}" ]] && [[ "${sample}" != "Input" ]]; then
		printf "${line}\t${tissue}\t${sample}\t${paired}\t${ref_dir}\n" >> $analysisfile
		sample_list+=("$shortname")
	fi
done < $samplefile

if [[ ${new_sample} != 0 ]]; then
#### Wait for the mapping sample scripts to finish
	printf "\nWaiting for the samples to be mapped...\n"
	wait ${pids[*]}

	i=0
	for check in ${check_list[@]}
	do
		if [ ! -e $check ]; then
			printf "\nProblem during mapping of ${checkname_list[i]}!\nCheck log: ${checkdatatype_list[i]}/logs/${checkname_list[i]}.log\n"
			exit 1
		else
			printf "\nSample ${checkname_list[i]} successfully mapped.\n"
		fi
		i=$((i+1))
	done
fi

### Plotting mapping stats for better visualization

printf "\nExtracting mapping stats table\n"
if [ -e combined/reports/temp_mapping_stats_${samplename}.txt ]; then
	rm -f combined/reports/temp_mapping_stats_${samplename}.txt
fi

while read line tissue sample rep sampleID path paired ref
do
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
	esac
	awk -v a=$line -v b=$tissue -v c=$sample -v d=$rep -v e=$ref '$1==a && $2==b && $3==c && $4==d && $5==e' ${datatype}/reports/summary_mapping_stats.txt >> combined/reports/temp_mapping_stats_${samplename}.txt
done < $samplefile

printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > combined/reports/summary_mapping_stats_${samplename}.txt
sort combined/reports/temp_mapping_stats_${samplename}.txt -u >> combined/reports/summary_mapping_stats_${samplename}.txt
rm -f combined/reports/temp_mapping_stats_${samplename}.txt
printf "\nPlotting mapping stats for all samples in the samplefile with R:\n"
R --version
Rscript --vanilla ${mc_dir}/MaizeCode_R_mapping_stats.r combined/reports/summary_mapping_stats_${samplename}.txt ${samplename}

if [[ "$keepgoing" == "STOP" ]]; then
	enddate=`date +%s`
	runtime=$((enddate-startdate))
	hours=$((runtime / 3600))
	minutes=$(( (runtime % 3600) / 60 ))
	seconds=$(( (runtime % 3600) % 60 ))
	printf "\nScript finished successfully without analysis in $hours:$minutes:$seconds (hh:mm:ss)!\n"		
	exit 0
fi	

#############################################################################################
########################################### PART3 ###########################################
############################# Starting analysis if not stopped ##############################
#############################################################################################

uniq_ref_list=($(printf "%s\n" "${ref_list[@]}" | sort -u))

if [ -e all_genes.txt ]; then
	rm -f all_genes.txt
fi

for ref in ${uniq_ref_list[@]}
do
	if [ -s ChIP/tracks/${ref}_all_genes.bed ]; then
		printf "ChIP/tracks/${ref}_all_genes.bed\n" >> all_genes.txt
	elif [ -s RNA/tracks/${ref}_all_genes.bed ]; then
		printf "RNA/tracks/${ref}_all_genes.bed\n" >> all_genes.txt
	else
		printf "Problem: no region file found for ${ref}!\n"
		exit 1
	fi
done

pids=()

if [[ "$wholeanalysis" == "STOP" ]]; then
	printf "\nPerforming only the single sample analysis\n"
	qsub -sync y -N maizecodeanalysis -o maizecode.log ${mc_dir}/MaizeCode_analysis.sh -f $analysisfile -r all_genes.txt -s &
	analysisname="${samplename}_no_region"
	check="combined/chkpts/${analysisname}"
else
	printf "\nPerforming the complete analysis on all genes\n"
	qsub -sync y -N maizecodeanalysis -o maizecode.log ${mc_dir}/MaizeCode_analysis.sh -f $analysisfile -r all_genes.txt &
	analysisname="${samplename}_on_all_genes"
	check="combined/chkpts/${analysisname}"
fi	
pids+=("$!")

wait ${pids[*]}

if [ -e all_genes.txt ]; then
	rm -f all_genes.txt
fi

if [ ! -e ${check} ]; then
	printf "\nProblem during the $analysisname analysis. Check maizecode.log\n"
	exit 1
else
	enddate=`date +%s`
	runtime=$((enddate-startdate))
	hours=$((runtime / 3600))
	minutes=$(( (runtime % 3600) / 60 ))
	seconds=$(( (runtime % 3600) % 60 ))
	printf "\nMaizeCode script finished successfully in $hours:$minutes:$seconds (hh:mm:ss)!\nCheck out the plots!\n"
fi

############################################################################################
########################################## MISC ############################################
############################################################################################

# # ##### To make the samplefile (e.g. for B73 endosperm)

# # printf "B73\tendosperm\tInput\tRep1\tS01\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me1\tRep1\tS02\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me3\tRep1\tS03\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K27ac\tRep1\tS04\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tInput\tRep2\tS05\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me1\tRep2\tS06\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K4me3\tRep2\tS07\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tH3K27ac\tRep2\tS08\t/seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846\tPE\tB73_v4\nB73\tendosperm\tRNAseq\tRep1\tendo\t/mnt/grid/martienssen/nlsas_norepl/data/data/archive/data/maizecode/released/B73/Long_Rampage/endo_rep1/RNAseq\tPE\tB73_v4\nB73\tendosperm\tRAMPAGE\tRep1\tendo\t/mnt/grid/martienssen/nlsas_norepl/data/data/archive/data/maizecode/released/B73/Long_Rampage/endo_rep1/RAMPAGE\tPE\tB73_v4\nB73\tendosperm\tRNAseq\tRep2\tendo\t/mnt/grid/martienssen/nlsas_norepl/data/data/archive/data/maizecode/released/B73/Long_Rampage/endo_rep2/RNAseq\tPE\tB73_v4\nB73\tendosperm\tRAMPAGE\tRep2\tendo\t/mnt/grid/martienssen/nlsas_norepl/data/data/archive/data/maizecode/released/B73/Long_Rampage/endo_rep2/RAMPAGE\tPE\tB73_v4\n" > B73_endosperm_samplefile.txt

############################################################################################
