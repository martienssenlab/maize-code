#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 8
#$ -l m_mem_free=12G
#$ -l tmp_free=100G
#$ -o ChIPanalysis.log
#$ -j y
#$ -N ChIPanalysis

usage="
##### Script for Maize code ChIP data analysis
#####
##### sh MaiCode_ChIP_analysis.sh -f samplefile [-h]
#####	-f: samplefile containing the samples to compare and in 4 tab-delimited columns:
##### 		Line, Tissue, Mark, PE or SE
##### 	-h: help, returns usage
##### 
##### It merges the two replicate files, and create pseudo-replicates by splitting the merged bam file into 2 halves
##### It calls broad (for H3K4me1) or narrow (for H3K4me3 and H3K27ac) peaks with Macs2 on each biological replicate,
##### on the merged file, and on both pseudo-replicates and will keep the peaks that are present in merged and both pseudo-replicates for downstream analysis,
##### it then creates bigwig files (log2 FC vs Input) for each biological rep and the merged file, makes fingerprint plots for each,
##### makes idr analysis between biological replicates, and calculates peak stats
#####
##### Requirements: samtools, bedtools, deeptools, macs2, idr
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

while getopts ":f:sh" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
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

tmp1=${samplefile##*temp_}
samplename=${tmp1%%_ChIP*}

if [ ! -s peaks/summary_peaks_${samplename}.txt ]; then
	printf "Line\tTissue\tMark\tPeaks_in_rep1\tPeaks_in_Rep2\tCommon_peaks\tCommon_peaks_IDR<=0.05\tPeaks_in_merged\tPeaks_in_pseudo_reps\tSelected_peaks\n" > peaks/summary_peaks_${samplename}.txt
fi

while read line tissue mark paired
do
	#### To merge bam files of replicates
	name=${line}_${tissue}_${mark}
	input=${line}_${tissue}_Input
	if [ ! -s mapped/${input}_merged.bam ]; then
		printf "\nMerging replicates of $input\n"
		samtools merge -@ $threads mapped/temp_${input}.bam mapped/rmdup_${input}_Rep1.bam mapped/rmdup_${input}_Rep2.bam
		samtools sort -@ $threads -o mapped/${input}_merged.bam mapped/temp_${input}.bam
		rm -f mapped/temp_${input}.bam
		samtools index -@ $threads mapped/${input}_merged.bam
	fi
	if [ ! -s mapped/${name}_merged.bam ]; then
		printf "\nMerging replicates of $name\n"
		samtools merge -f -@ $threads mapped/temp_${name}.bam mapped/rmdup_${name}_Rep1.bam mapped/rmdup_${name}_Rep2.bam
		samtools sort -@ $threads -o mapped/${name}_merged.bam mapped/temp_${name}.bam
		rm -f mapped/temp_${name}.bam
		samtools index -@ $threads mapped/${name}_merged.bam
	fi
	if [ ! -s mapped/${name}_pseudo1.bam ]; then
		printf "\nSplitting $name in two pseudo-replicates\n"
		samtools view -b -h -s 1.5 -@ $threads -U mapped/temp_${name}_pseudo2.bam -o mapped/temp_${name}_pseudo1.bam mapped/${name}_merged.bam
		samtools sort -@ $threads -o mapped/${name}_pseudo1.bam mapped/temp_${name}_pseudo1.bam
		rm -f mapped/temp_${name}_pseudo1.bam
		samtools index -@ $threads mapped/${name}_pseudo1.bam
		samtools sort -@ $threads -o mapped/${name}_pseudo2.bam mapped/temp_${name}_pseudo2.bam
		rm -f mapped/temp_${name}_pseudo2.bam
		samtools index -@ $threads mapped/${name}_pseudo2.bam
	fi
	#### To call either broad or narrow peaks if not already exisiting
	case "$mark" in
		H3K4me1) peaktype="broad";;
		H3K4me3) peaktype="narrow";;
		H3K27ac) peaktype="narrow";;
	esac
	pidsb=()
	for filetype in merged Rep1 Rep2 pseudo1 pseudo2
	do
		case "$filetype" in
			Rep1|Rep2) 	namefiletype=mapped/rmdup_${name}_${filetype}.bam
						inputfiletype=mapped/rmdup_${input}_${filetype}.bam
						param="";;
			pseudo1|pseudo2)	namefiletype=mapped/${name}_${filetype}.bam
								inputfiletype=mapped/${input}_merged.bam
								param="";;
			merged)	namefiletype=mapped/${name}_${filetype}.bam
					inputfiletype=mapped/${input}_${filetype}.bam
					param="-B";;
		esac		
		printf "\nStarting single ChIP sample analysis for $name $filetype\n"
		qsub -N ${name}_${filetype} -V -cwd -sync y -pe threads 20 -l m_mem_free=8G -l tmp_free=50G -j y -o logs/analysis_${name}_${filetype}.log <<-EOF &
			#!/bin/bash
			set -e -o pipefail
		
			export threads=$NSLOTS
		
			if [[ $paired == "PE" ]]; then
				if [[ $peaktype == "broad" ]] && [ ! -s peaks/${name}_${filetype}_peaks.${peaktype}Peak ]; then
					printf "\nCalling broad peaks for PE $namefiletype (vs $inputfiletype) with macs2 version:\n"
					macs2 --version
					macs2 callpeak -t ${namefiletype} -c ${inputfiletype} -f BAMPE -g 2.2e9 ${param} -n ${name}_${filetype} --keep-dup "all" --outdir peaks/ --tempdir $TMPDIR --broad
				elif [[ $peaktype == "narrow" ]] && [ ! -s peaks/${namefiletype}_peaks.${peaktype}Peak ]; then
					printf "\nCalling narrow peaks for PE $namefiletype (vs $inputfiletype) with macs2 version:\n"
					macs2 --version
					macs2 callpeak -t ${namefiletype} -c ${inputfiletype} -f BAMPE -g 2.2e9 ${param} -n ${name}_${filetype} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR
				elif [ -s peaks/${name}_${filetype}_peaks.${peaktype}Peak ]; then
					printf "\n$peaktype peaks already called for $namefiletype\n"
				else
					printf "\nSomething is wrong! Check usage by running the script without arguments\n"
					printf "$usage\n"
					exit 1
				fi
			elif [[ $paired == "SE" ]]; then
				if [[ $peaktype == "broad" ]] && [ ! -s peaks/${name}_${filetype}_peaks.${peaktype}Peak ]; then
					printf "\nCalling broad peaks for SE $namefiletype (vs $inputfiletype) with macs2 version:\n"
					macs2 --version
					macs2 callpeak -t ${namefiletype} -c ${inputfiletype} -f BAM -g 2.2e9 ${param} -n ${name}_${filetype} --keep-dup "all" --outdir peaks/ --tempdir $TMPDIR --nomodel --broad
				elif [[ $peaktype == "narrow" ]] && [ ! -s peaks/${name}_${filetype}_peaks.${peaktype}Peak ]; then
					printf "\nCalling narrow peaks for SE $namefiletype (vs $inputfiletype) with macs2 version:\n"
					macs2 --version
					macs2 callpeak -t ${namefiletype} -c ${inputfiletype} -f BAM -g 2.2e9 ${param} -n ${name}_${filetype} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR --nomodel
				elif [ -s peaks/${name}_${filetype}_peaks.${peaktype}Peak ]; then
					printf "\n$peaktype peaks already called for $namefiletype\n"
				else
					printf "\nSomething is wrong! Check usage:\n"
					printf "$usage\n"
					exit 1
				fi
			else
				printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
				exit 1
			fi
			#### To create bw files if not already exisiting
			if [ ! -s tracks/${name}_${filetype}.bw ]; then
				printf "\nMaking bigwig files for $namefiletype with deeptools version:\n"
				deeptools --version
				bamCompare -b1 ${namefiletype} -b2 ${inputfiletype} -o tracks/${name}_${filetype}.bw -p $threads --binSize 1 --scaleFactorsMethod "None" --normalizeUsing CPM
			else
				printf "\nBigwig file for $namefiletype already exists\n"
			fi
			#### To create fingerprint plots if not already exisiting
			if [ ! -s plots/Fingerprint_${name}_${filetype}.png ]; then
				printf "\nPlotting Fingerprint for $namefiletype with deeptools version:\n"
				deeptools --version
				plotFingerprint -b ${namefiletype} ${inputfiletype} -o plots/Fingerprint_${name}_${filetype}.png -p $threads -l ${name} ${input}
			else
				printf "\nFringerprint plot for $namefiletype already exists\n"
			fi
		EOF
		pidsb+=("$!")
	done
	printf "\nWaiting for $name single sample files to be processed\n"
	wait ${pidsb[*]}
	
	#### To get IDR analysis on biological replicates
	if [ ! -s peaks/idr_${name}.${peaktype}Peak ]; then
		printf "\nDoing IDR analysis on both replicates from ${line}_${tissue}_${mark} ($peaktype peaks) with idr version:\n"
		idr --version
		idr --input-file-type ${peaktype}Peak --output-file-type ${peaktype}Peak --samples peaks/${name}_Rep1_peaks.${peaktype}Peak peaks/${name}_Rep2_peaks.${peaktype}Peak -o peaks/idr_${name}.${peaktype}Peak -l reports/idr_${name}.log --plot
	else
		printf "\nIDR analysis already done for ${name}\n"
	fi
	#### To get the final selected peak file (peaks called in merged also present in both replicates)
	awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_merged.${peaktype}Peak | sort -k1,1 -k2,2n -u > peaks/temp_${name}.bed
	bedtools intersect -a peaks/${name}_pseudo1.bed -b peaks/${name}_pseudo2.bed > peaks/temp2_${name}.bed
	bedtools intersect -a peaks/temp_${name}.bed -b peaks/temp2_${name}.bed -u > peaks/temp3_${name}.bed
	bedtools intersect -a peaks/${name}_merged.${peaktype}Peak -b peaks/temp3_${name}.bed -u > peaks/selected_peaks_${name}.${peaktype}Peak

	#### To get some peaks stats for each mark
	printf "\nCalculating peak stats for ${name} in ${peaktype} peaks\n"
	rep1=$(awk '{print $1,$2,$3}' peaks/${name}_Rep1_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
	rep2=$(awk '{print $1,$2,$3}' peaks/${name}_Rep2_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
	common=$(awk '{print $1,$2,$3}' peaks/idr_${name}.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
	idr=$(awk '$5>=540 {print $1,$2,$3}' peaks/idr_${name}.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
	merged=$(awk '{print $1,$2,$3}' peaks/${name}_merged_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
	pseudos=$(awk '{print $1,$2,$3}' peaks/temp2_${name}.bed | sort -k1,1 -k2,2n -u | wc -l)
	selected=$(awk '{print $1,$2,$3}' peaks/selected_peaks_${name}.bed | sort -k1,1 -k2,2n -u | wc -l)
	awk -v OFS="\t" -v a=$line -v b=$tissue -v c=$mark -v d=$rep1 -v e=$rep2 -v f=$common -v g=$idr -v h=$merged -v i=$pseudos -v j=$selected 'BEGIN {print a,b,c,d,e,f" ("f/d*100"%rep1;"f/e*100"%rep2)",g" ("g/f*100"%common)",h,i,j" ("j/h*100"%merged)"}' >> peaks/summary_peaks_${samplename}.txt
	rm -f peaks/temp*
	touch chkpts/analysis_${name}
done < $samplefile

printf "\nScript finished successfully!\n"

