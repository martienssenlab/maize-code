#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=1.5G
#$ -l tmp_free=20G
#$ -o logs/mCsample.log
#$ -j y
#$ -N mCsample

usage="
##### Script for Maize code DNA methylation data analysis, used by script MaizeCode.sh for mC samples
#####
##### sh MaizeCode_mC_sample.sh -x datatype -d reference directory -l inbred line -t tissue -m histone mark -r replicate ID -i sample ID -f path to sample -p paired -s step -a mappingoption
##### 	-x: type of data (should be 'mC')
##### 	-d: folder containing the reference directory (e.g. ~/data/Genomes/Zea_mays/B73_v4)
##### 	-l: inbred line (e.g. B73)
##### 	-t: tissue (e.g. endosperm)
##### 	-m: Library method [ Pico | EMseq | mC ] If it is not one of these options, it defaults to mC
##### 	-r: replicate ID (e.g. Rep1)
#####	-i: sample ID (name in original folder or SRR number)
#####	-f: path to original folder or SRA
##### 	-p: if data is paired-end (PE) or single-end (SE) [ PE | SE ]
#####	-s: status of the raw data [ download | trim | done ] 'download' if sample needs to be copied/downloaded, 'trim' if only trimming has to be performed, 'done' if trimming has already been performed
#####	-a: what option to use for mapping [ default | all | colcen | colcenall | heavy ] (colcen: very-sensitive, -k 100; all: no MAPQ>10; heavy: RNA only, default here)
##### 	-h: help, returns usage
#####
##### It downloads or copies the files, runs fastQC, trims adapters with cutadapt, aligns, deduplicate and extract methylation with bismark and get some mapping stats
##### Pico changes the directionality of the libraries. If lambda was found in the reference DNA sequences it is used for non-conversion rate, otherwise plastid DNA will be used if found.
#####
##### Requirements: samtools, fastQC, Cutadapt, Bowtie2, Bismark, parallel-fastq-dump (if downloading SRA data)
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS
export limthreads=$((threads/3))

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
		m)	export met=${OPTARG};;
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

if [ ! ${data} ] || [ ! ${ref_dir} ] || [ ! ${line} ] || [ ! ${tissue} ] || [ ! ${met} ] || [ ! ${rep} ] || [ ! ${sampleID} ] || [ ! ${path} ] || [ ! ${paired} ] || [ ! ${step} ]; then
	printf "Missing arguments!\n"
	printf "${usage}\n"
	exit 1
fi

if [ ! ${mapparam} ]; then
	printf "No mapping option selected, using default\n"
	export mapparam="default"
elif [[ "${mapparam}" == "default" || "${mapparam}" == "colcen" || "${mapparam}" == "colcenall" || "${mapparam}" == "all" || "${mapparam}" == "heavy" ]]; then
	printf "${mapparam} chosen as the mapping option\n"
else
	printf "Unknown mapping option selected\n"
	printf "${usage}\n"
	exit 1
fi

export ref=${ref_dir##*/}

name=${line}_${tissue}_mC_${rep}

if [[ ${paired} == "PE" ]]; then
	if [[ ${met} == "Pico" ]]; then
		param1="-u 10 -U 10 -q 10 -m 20"
		param2="--non_directional --maxins 1000"
		param3=""
	elif [[ ${met} == "EMseq" ]]; then
		param1="-u 10 -U 10 -q 10 -m 20"
		param2="--maxins 1000"
		param3=""
	else
		param1="-q 10 -m 20"
		param2="--maxins 1000" 
		param3="--ignore_r2 2"
	fi
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
   		cutadapt -j ${threads} ${param1} -a AGATCGGAAGAGCACACGTCTGAAC -A AGATCGGAAGAGCGTCGTGTAGGGA -o fastq/trimmed_${name}_R1.fastq.gz -p fastq/trimmed_${name}_R2.fastq.gz fastq/${name}_R1.fastq.gz fastq/${name}_R2.fastq.gz |& tee reports/trimming_${name}.txt
		#### Removing untrimmed fastq
		rm -f fastq/${name}_R*.fastq.gz
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for ${name}\n"
		fastqc -o reports/ fastq/trimmed_${name}_R1.fastq.gz
		fastqc -o reports/ fastq/trimmed_${name}_R2.fastq.gz
	fi
	#### Proceeding with DNA methylation analysis with Bismark
	R1="fastq/trimmed_${name}_R1.fastq.gz"
	R2="fastq/trimmed_${name}_R2.fastq.gz"
	shortR1="trimmed_${name}_R1"
	printf "\nAligning ${name} with bismark/bowtie2\n"
	bismark --genome ${ref_dir} ${param2} --local --multicore ${limthreads} --temp_dir=${TMPDIR} -o mapped/${name} --gzip --nucleotide_coverage -1 ${R1} -2 ${R2} |& tee reports/alignment_bismark_${name}.txt
	printf "\nDeduplicating with bismark\n"
	deduplicate_bismark -p --output_dir mapped/${name}/ -o ${name} --bam mapped/${name}/${shortR1}_bismark_bt2_pe.bam |& tee reports/deduplication_bismark_${name}.txt
	rm -f mapped/${name}/${shortR1}_bismark_bt2_pe.bam
	printf "\nCalling mC for ${name}"
	bismark_methylation_extractor -p --comprehensive -o methylcall/ ${param3} --gzip --multicore ${limthreads} --buffer_size 10G --cytosine_report --CX --genome_folder ${ref_dir} mapped/${name}/${name}.deduplicated.bam
	rm -f methylcall/C*context_${name}*
	rm -f methylcall/${name}*bismark.cov*
	printf "\nMaking final html report for ${name}\n"
	bismark2report -o final_report_${name}.html --dir reports/ --alignment_report mapped/${name}/trimmed_${name}_R1_bismark_bt2_PE_report.txt --dedup_report mapped/${name}/trimmed_${name}_R1_bismark_bt2_pe.deduplication_report.txt --splitting_report methylcall/${name}.deduplicated_splitting_report.txt --mbias_report methylcall/${name}.deduplicated.M-bias.txt --nucleotide_report mapped/${name}/trimmed_${name}_R1_bismark_bt2_pe.nucleotide_stats.txt
 	printf "\nCalculting coverage stats for ${name}\n"
	tot=$(cat reports/alignment_bismark_${name}.txt | grep "Sequence pairs analysed in total:" | awk -v FS="\t" 'END {print $2}')
	map=$(cat reports/alignment_bismark_${name}.txt | grep "Number of paired-end alignments with a unique best hit:" | awk -v FS="\t" 'END {print $2}')
  	uniq=$(cat reports/deduplication_bismark_${name}.txt | grep "Total count of deduplicated leftover sequences:" | awk -v FS=":" 'END {print $2}' | awk '{print $1}')
  	if grep -E -q "J02459.1_48502" ${ref_dir}/chrom.sizes; then
    		zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v l=${line} -v t=${tissue} -v r=${rep} -v z=${tot} -v y=${map} -v x=${uniq} '{a+=1; b=$4+$5; g+=b; if ($1=="J02459.1_48502") {m+=$4; n+=b;}; if (b>0) {c+=1; d+=b;} else f+=1; if (b>2) e+=1} END {print l,t,r,z,y,x,c/a*100,e/a*100,g/a,d/c,m/n*100}' >> reports/summary_mapping_stats.txt
  	elif grep -E -q "Pt|ChrC|chrc" ${ref_dir}/chrom.sizes; then
  		zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v l=${line} -v t=${tissue} -v r=${rep} -v z=${tot} -v y=${map} -v x=${uniq} '{a+=1; b=$4+$5; g+=b; if ($1 == "Pt" || $1 == "ChrC" || $1 == "chrC") {m+=$4; n+=b;}; if (b>0) {c+=1; d+=b;} else f+=1; if (b>2) e+=1} END {print l,t,r,z,y,x,c/a*100,e/a*100,g/a,d/c,m/n*100}' >> reports/summary_mapping_stats.txt
  	else
    		zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v l=${line} -v t=${tissue} -v r=${rep} -v z=${tot} -v y=${map} -v x=${uniq} '{a+=1; b=$4+$5; g+=b; if (b>0) {c+=1; d+=b;} else f+=1; if (b>2) e+=1} END {print l,t,r,z,y,x,c/a*100,e/a*100,g/a,d/c,"NA"}' >> reports/summary_mapping_stats.txt
  	fi
	zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v s=${name} '($4+$5)>0 {a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHH.bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHG.bedGraph"; else print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CG.bedGraph"}'
	for strand in plus minus
	do
		case "${strand}" in 
			plus)	sign="+";;
			minus)	sign="-";;
		esac
		zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v n=${sign} '$3==n' | awk -v OFS="\t" -v s=${name} -v d=${strand} '($4+$5)>0 {a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHH_"d".bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHG_"d".bedGraph"; else if ($6=="CG") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CG_"d".bedGraph"}'
   	done
	for context in CG CHG CHH
	do
		printf "\nMaking bigwig files of ${context} context for ${name}\n"
		LC_COLLATE=C sort -k1,1 -k2,2n methylcall/${name}_${context}.bedGraph > methylcall/sorted_${name}_${context}.bedGraph
		bedGraphToBigWig methylcall/sorted_${name}_${context}.bedGraph ${ref_dir}/chrom.sizes methylcall/${name}_${context}.bw
   		for strand in plus minus
      		do
 			printf "\nMaking ${strand} strand bigwig files of ${context} context for ${name}\n"
			LC_COLLATE=C sort -k1,1 -k2,2n methylcall/${name}_${context}_${strand}.bedGraph > methylcall/sorted_${name}_${context}_${strand}.bedGraph
			bedGraphToBigWig methylcall/sorted_${name}_${context}_${strand}.bedGraph ${ref_dir}/chrom.sizes methylcall/${name}_${context}_${strand}.bw
    		done
	done
	rm -f methylcall/*${name}*bedGraph*
elif [[ ${paired} == "SE" ]]; then
	if [[ ${met} == "Pico" ]]; then
		param1="-u 10 -q 10 -m 20"
		param2="--non_directional"
		param3=""
	elif [[ ${met} == "EMseq" ]]; then
		param1="-u 10 -q 10 -m 20"
		param2=""
		param3=""
	else
		param1="-q 10 -m 20"
		param2="" 
		param3=""
	fi
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
		cutadapt -j ${threads} ${param1} -a AGATCGGAAGAGCACACGTCTGAAC -o fastq/trimmed_${name}.fastq.gz fastq/${name}.fastq.gz |& tee reports/trimming_${name}.txt
		#### Removing untrimmed fastq
		rm -f fastq/${name}.fastq.gz
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for ${name}\n"
		fastqc -o reports/ fastq/trimmed_${name}.fastq.gz
	fi
	#### Proceeding with DNA methylation analysis with Bismark
	printf "\nAligning ${name} with bismark_bowtie2\n"
	bismark --genome ${ref_dir} ${param2} --local --multicore ${limthreads} --temp_dir=${TMPDIR} -o mapped/${name} --gzip --nucleotide_coverage fastq/trimmed_${name}.fastq.gz |& tee reports/alignment_bismark_${name}.txt
	printf "\nDeduplicating ${name} with bismark\n"
	deduplicate_bismark -s --output_dir mapped/${name}/ -o ${name} --bam mapped/${name}/trimmed_${name}_bismark_bt2.bam |& tee reports/deduplication_bismark_${name}.txt
	printf "\nCalling mC for ${name}"
	bismark_methylation_extractor -s --comprehensive -o methylcall/ ${param3} --gzip --multicore ${limthreads} --buffer_size 10G --cytosine_report --CX --genome_folder ${ref_dir} mapped/${name}/${name}.deduplicated.bam
	rm -f methylcall/C*context_${name}*
	rm -f methylcall/${name}*bismark.cov*
	printf "\nMaking final html report for ${name}\n"
	bismark2report -o final_report_${name}.html --dir reports/ --alignment_report mapped/${name}/trimmed_${name}_bismark_bt2_SE_report.txt --dedup_report mapped/${name}/trimmed_${name}_bismark_bt2.deduplication_report.txt --splitting_report methylcall/${name}.deduplicated_splitting_report.txt --mbias_report methylcall/${name}.deduplicated.M-bias.txt --nucleotide_report mapped/${name}/trimmed_${name}_bismark_bt2.nucleotide_stats.txt
 	# printf "\nCalculting coverage stats for ${name}\n" To be completed	
	zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v s=${name} '($4+$5)>0 {a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHH.bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHG.bedGraph"; else print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CG.bedGraph"}'
	for strand in plus minus
	do
		case "${strand}" in 
			plus)	sign="+";;
			minus)	sign="-";;
		esac
		zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v n=${sign} '$3==n' | awk -v OFS="\t" -v s=${name} -v d=${strand} '($4+$5)>0 {a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHH_"d".bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHG_"d".bedGraph"; else if ($6=="CG") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CG_"d".bedGraph"}'
   	done
	for context in CG CHG CHH
	do
		printf "\nMaking bigwig files of ${context} context for ${name}\n"
		LC_COLLATE=C sort -k1,1 -k2,2n methylcall/${name}_${context}.bedGraph > methylcall/sorted_${name}_${context}.bedGraph
		bedGraphToBigWig methylcall/sorted_${name}_${context}.bedGraph ${ref_dir}/chrom.sizes methylcall/${name}_${context}.bw
   		for strand in plus minus
      		do
 			printf "\nMaking ${strand} strand bigwig files of ${context} context for ${name}\n"
			LC_COLLATE=C sort -k1,1 -k2,2n methylcall/${name}_${context}_${strand}.bedGraph > methylcall/sorted_${name}_${context}_${strand}.bedGraph
			bedGraphToBigWig methylcall/sorted_${name}_${context}_${strand}.bedGraph ${ref_dir}/chrom.sizes methylcall/${name}_${context}_${strand}.bw
    		done
	done
	rm -f methylcall/*${name}*bedGraph*
else
	printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
	exit 1
fi

printf "\nScript finished successfully!\n"
touch chkpts/${name}_${ref}
