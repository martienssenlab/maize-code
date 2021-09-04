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
##### sh MaizeCode.sh -f <samplefile> -p <path to genome reference> [-s] [-c] [-t] [-h]
##### 	-f: samplefile
##### 	-p: path to the folder containing all the different genome references (e.g. ~/data/Genomes/Zea_mays)
#####	-s: if set, the whole analysis does not proceed (default=not set, keep going with the analysis over all the samples in the samplefile)
#####	-c: if set, only single samples analysis proceeds, not grouped analysis per line (default=not set, keep going with the complete analysis)
#####	-t: if set, only partial grouped analysis per line, no heatmaps with deeptools (default=not set, keep going with the complete analysis)
##### 	-h: help, returns usage
#####
##### The samplefile should be a tab-delimited text file with 8 columns:
##### col #1: Type of data [ RNAseq | RAMPAGE | ChIP | TF_*] (shRNA in development). All options are case-sensitive. For TF_*, the star should be replaced by the name of the TF, e.g. TF_TB1.
##### col #2: Line (e.g. B73)
##### col #3: Tissue (e.g endosperm) 
##### col #4: Sample (e.g. 'H3K4me3' or 'Input' for ChIP, 'IP' or 'Input' for TF CHIPseq, shRNA, RNAseq or RAMPAGE for RNA (same as data type).
##### col #5: Replicate ID [ Rep1 | Rep2 ] (can be more for RNA samples, not for ChIP/TF where it can only be 1 or 2).
##### col #6: SequencingID (e.g. S01). Unique identifier for the name of the sample in the raw sequencing folder which path is given in the next column. If downloading from SRA, put the SRR ID here.
##### col #7: Path to the fastq files (e.g. /seq/Illumina_runs/NextSeqData/NextSeqOutput/190913_NB501555_0636_AH5HG7BGXC/Data/Intensities/BaseCalls/304846). If downloading from SRA, put 'SRA'.
##### col #8: If data is paired-end or single-end [ PE | SE ]. 
##### col #9: Name of the genome reference to map (e.g. B73_v4). Each genome reference should have a unique folder that contains a single fasta file, gff3 file and gtf file (can all be gzipped).
##### The gff3 files should have 'gene' in column 3 and exons should be linked by 'Parent' in column 9
##### The fasta and gff3 files should have the same chromosome names (i.e. 1 2 3... and 1 2 3... or Chr1 Chr2 Chr3... and Chr1 Chr2 Chr3...)
##### For cleaner naming purposes, use '_samplefile.txt' as suffix
#####
##### This script creates the folders needed,
##### prepares the genome index and other environment files if they are not done for each type of data,
##### sends each sample into a datatype-specific script for mapping (MaizeCode_${datatype}_sample.sh) if they have not been mapped before
##### starts the script (MaizeCode_analysis.sh) for the analysis on all the samples in the samplefile (if -s is not set).
#####
"

set -e -o pipefail

printf "\n\n"
date
startdate=`date +%s`
printf "\n"

export threads=${NSLOTS}
export mc_dir="${PWD}/scripts/"
printf "\nRunning MaizeCode.sh script in working directory ${PWD}\n"

if [ $# -eq 0 ]; then
	printf "${usage}\n"
	exit 1
fi

while getopts "f:p:scth" opt; do
	case $opt in
		h) 	printf "${usage}\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		p)	export pathtoref=${OPTARG};;
		s)	printf "\nOption not to perform analysis selected\n"
			export keepgoing="STOP";;
		c)	printf "\nOption not to perform combined analysis selected\n"
			export wholeanalysis="STOP";;
		t)	printf "\nOption to perform partial combined analysis selected\n"
			export total="No";;
		*)	printf "\nUnknown option\n${usage}\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! ${samplefile} ]; then
	printf "Samplefile missing!\n"
	printf "${usage}\n"
	exit 1
fi
if [ ! ${pathtoref} ]; then
	printf "Path to reference genome folders missing!\n"
	printf "${usage}\n"
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
while read data line tissue sample rep sampleID path paired ref
do
	name=${line}_${tissue}_${sample}_${rep}
	case "${data}" in
		ChIP) 	env="ChIP";;
		RNAseq) env="RNA";;
		RAMPAGE) env="RNA";;
		shRNA) env="shRNA";;
		TF_*) env="TF";;
		*) env="unknown";;
	esac
	if [[ "${env}" == "unknown" ]]; then
		printf "Type of data unknown!\n${usage}\n"
		exit 1
	fi
	if [ ! -e ${env}/chkpts/env_${ref} ]; then
		new_env=1
		newdatatype_list+=("${env}")
		newref_list+=("${ref}")
		data_ref_list+=("${env}_${ref}")
	fi
done < ${samplefile}

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
		for env in ${uniq_newdatatype_list[@]}
		do
			if [[ " ${data_ref_list[@]} " =~ " ${env}_${ref} " ]]; then
				check_list+=("${env}/chkpts/env_${ref}")
				if [ ! -d ./${env} ]; then
					mkdir ./${env}
					mkdir ./${env}/fastq
					mkdir ./${env}/mapped
					mkdir ./${env}/tracks
					mkdir ./${env}/reports
					mkdir ./${env}/logs
					mkdir ./${env}/chkpts
					mkdir ./${env}/plots
				fi
				printf "\nPreparing environment of ${ref} genome for ${env} data\n"
				qsub -sync y -N env_${ref}_${env} -o ${env}/logs/env_${ref}.log ${mc_dir}/MaizeCode_check_environment.sh -p ${pathtoref} -r ${ref} -d ${env} &
				pids+=("$!")
			elif [[ ! " ${data_ref_list[@]} " =~ " ${env}_${ref} " ]]; then
			### Combination folder * datatype * ref does not exist in the sample file, moving on
				:
			elif [ ! -d ${pathtoref}/${ref} ]; then
				printf "\nNo ${ref} folder found in ${pathtoref}\n"
				printf "${usage}\n"
				exit 1
			else
				printf "\nError!\n"
				printf "${usage}\n"
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
		env=${check%%/*}
		ref=${check##*/env_}
		if [ ! -e ${check} ]; then
			printf "\nProblem while making ${env} environment with ${ref} genome!\nCheck log: ${env}/logs/env_${ref}.log\n"
			exit 1
		else
			printf "\nEnvironment for ${env} with ${ref} genome now good to go\n"
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

new_sample=0
check_list=()
checkname_list=()
checkdatatype_list=()
ref_list=()
sample_list=()
pids=()
while read data line tissue sample rep sampleID path paired ref
do
	ref_dir=${pathtoref}/${ref}
	case "${data}" in
		ChIP) 	env="ChIP"
			shortname=${line}_${tissue}_${sample}
			name=${line}_${tissue}_${sample}_${rep};;
		RNAseq) env="RNA"
			shortname=${line}_${tissue}_${sample}
			name=${line}_${tissue}_${sample}_${rep};;
		RAMPAGE) 	env="RNA"
				shortname=${line}_${tissue}_${sample}
				name=${line}_${tissue}_${sample}_${rep};;
		shRNA) 	env="shRNA"
			shortname=${line}_${tissue}_${sample}
			name=${line}_${tissue}_${sample}_${rep};;
		TF_*) 	env="TF"
			tmp=${data##TF_}
			shortname=${line}_${tmp}_${sample}
			name=${line}_${tmp}_${sample}_${rep};;
	esac
	check=${env}/chkpts/${name}_${ref}
	ref_list+=("${ref}")
	if [ -e ${check} ]; then
		printf "Sample ${name} has already been mapped to ${ref} genome\n"
	else
		new_sample=1
		checkname_list+=("${name}")
		checkdatatype_list+=("${env}")
		check_list+=("${check}")
		if ls ./${env}/fastq/trimmed_${name}*.fastq.gz 1> /dev/null 2>&1; then
			printf "\nTrimmed fastq file(s) for ${name} already exist\n"
			export step="done"
		elif ls ./${env}/fastq/${name}*.fastq.gz 1> /dev/null 2>&1; then
			printf "\nFastq file(s) for ${name} already exist\n"
			export step="trim"
		else
			printf "\nNew sample ${name} to be copied/downloaded\n"
			export step="download"
		fi
		printf "\nRunning ${env} mapping script for ${name} on ${ref} genome\n"
		cd ${env}
		qsub -sync y -N ${name} -o logs/${name}.log ${mc_dir}/MaizeCode_${env}_sample.sh -x ${data} -d ${ref_dir} -l ${line} -t ${tissue} -m ${sample} -r ${rep} -i ${sampleID} -f ${path} -p ${paired} -s ${step} &
		pids+=("$!")
		cd ..
	fi
	if [[ ! "${sample_list[@]}" =~ "${shortname}" ]] && [[ "${sample}" != "Input" ]]; then
		printf "${data}\t${line}\t${tissue}\t${sample}\t${paired}\t${ref_dir}\n" >> ${analysisfile}
		sample_list+=("${shortname}")
	fi
done < ${samplefile}

if [[ ${new_sample} != 0 ]]; then
#### Wait for the mapping sample scripts to finish
	printf "\nWaiting for the samples to be mapped...\n"
	wait ${pids[*]}

	i=0
	for check in ${check_list[@]}
	do
		if [ ! -e ${check} ]; then
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
if [ -e combined/reports/temp2_mapping_stats_${samplename}.txt ]; then
	rm -f combined/reports/temp2_mapping_stats_${samplename}.txt
fi

while read data line tissue sample rep sampleID path paired ref
do
	case "${data}" in
		ChIP) env="ChIP"
			stat="plot1"
			name="${tissue}";;
		RNAseq) env="RNA"
			stat="plot1"
			name="${tissue}";;
		RAMPAGE) env="RNA"
			stat="plot1"
			name="${tissue}";;
		TF_*) env="TF"
			stat="plot1"
			name=${data##TF_};;
		shRNA) env="shRNA"
			stat="plot2"
			name="${line}_${tissue}_${sample}_${rep}";;
	esac
	if [[ ${stat} == "plot1" ]]; then
		awk -v a=${line} -v b=${name} -v c=${sample} -v d=${rep} -v e=${ref} '$1==a && $2==b && $3==c && $4==d && $5==e' ${env}/reports/summary_mapping_stats.txt >> combined/reports/temp_mapping_stats_${samplename}.txt
	elif [[ ${stat} == "plot2" ]]; then
		cat shRNA/reports/sizes_*${sample}* >> combined/reports/temp2_mapping_stats_${samplename}.txt
	fi
done < ${samplefile}

if [ -s combined/reports/temp_mapping_stats_${samplename}.txt ]; then
	printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > combined/reports/summary_mapping_stats_${samplename}_ChIP_RNA.txt
	sort combined/reports/temp_mapping_stats_${samplename}.txt -u >> combined/reports/summary_mapping_stats_${samplename}_ChIP_RNA.txt
	rm -f combined/reports/temp_mapping_stats_${samplename}.txt
	printf "\nPlotting mapping stats for all ChIP and RNA samples in the samplefile with R:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_mapping_stats.r combined/reports/summary_mapping_stats_${samplename}_ChIP_RNA.txt ${samplename}
fi
if [ -s combined/reports/temp2_mapping_stats_${samplename}.txt ]; then
	printf "Sample\tType\tSize\tCount\n" > combined/reports/summary_mapping_stats_${samplename}_shRNAs.txt
	sort combined/reports/temp2_mapping_stats_${samplename}.txt -u >> combined/reports/summary_mapping_stats_${samplename}_shRNA.txt
	rm -f combined/reports/temp2_mapping_stats_${samplename}.txt
	printf "\nPlotting mapping stats for all shRNA samples in the samplefile with R:\n"
	R --version
	Rscript --vanilla ${mc_dir}/MaizeCode_R_shRNA_stats.r combined/reports/summary_mapping_stats_${samplename}_shRNA.txt ${samplename}
fi 

if [[ "$keepgoing" == "STOP" ]]; then

	printf "\nPerforming multiQC analysis on all samples from ${samplename}\n"
	multiqc --version
	multiqc . -n multiqc_${samplename} -o combined/logs/ -i "multiQC_${samplename}" -f -ip -e deeptools -e snippy -z

	enddate=`date +%s`
	runtime=$((enddate-startdate))
	hours=$((runtime / 3600))
	minutes=$(( (runtime % 3600) / 60 ))
	seconds=$(( (runtime % 3600) % 60 ))
	printf "\nScript finished successfully without analysis in ${hours}:${minutes}:${seconds} (hh:mm:ss)!\n"		
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
	elif [ -s TF/tracks/${ref}_all_genes.bed ]; then
		printf "TF/tracks/${ref}_all_genes.bed\n" >> all_genes.txt
	elif [ -s shRNA/tracks/${ref}_all_genes.bed ]; then
		printf "shRNA/tracks/${ref}_all_genes.bed\n" >> all_genes.txt
	else
		printf "Problem: no region file found for ${ref}!\n"
		exit 1
	fi
done

pids=()

if [[ "$wholeanalysis" == "STOP" ]]; then
	printf "\nPerforming only the single sample analysis\n"
	qsub -sync y -N maizecodeanalysis -o maizecode.log ${mc_dir}/MaizeCode_analysis.sh -f ${analysisfile} -r all_genes.txt -s &
	analysisname="${samplename}_no_region"
	check="combined/chkpts/${analysisname}"
elif [[ "$total" == "No" ]]; then
	printf "\nPerforming partial analysis on all genes\n"
	qsub -sync y -N maizecodeanalysis -o maizecode.log ${mc_dir}/MaizeCode_analysis.sh -f ${analysisfile} -r all_genes.txt -t &
	analysisname="${samplename}_on_all_genes"
	check="combined/chkpts/${analysisname}"
else
	printf "\nPerforming complete analysis on all genes\n"
	qsub -sync y -N maizecodeanalysis -o maizecode.log ${mc_dir}/MaizeCode_analysis.sh -f ${analysisfile} -r all_genes.txt &
	analysisname="${samplename}_on_all_genes"
	check="combined/chkpts/${analysisname}"
fi	
pids+=("$!")

wait ${pids[*]}

printf "\nPerforming multiQC analysis on all samples from ${samplename}\n"
multiqc --version
multiqc . -n multiqc_${samplename} -o combined/logs/ -i "multiQC_${samplename}" -f -ip -e deeptools -e snippy -z

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
	printf "\nMaizeCode script finished successfully in ${hours}:${minutes}:${seconds} (hh:mm:ss)!\nCheck out the plots!\n"
fi
