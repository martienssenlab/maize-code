#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=3G
#$ -l tmp_free=20G
#$ -o ChIPanalysis.log
#$ -j y
#$ -N ChIPanalysis

usage="
##### Script for Maize code Histone ChIP data analysis, used by script MaizeCode_analysis.sh for ChIP data
#####
##### sh MaiCode_ChIP_analysis.sh -f samplefile [-h]
#####	-f: samplefile containing the samples to compare and the reference directory in 6 tab-delimited columns:
##### 		Data, Line, Tissue, Mark, PE or SE, Reference directory
##### 	-h: help, returns usage
##### 
##### It merges the two replicate files, and creates pseudo-replicates by splitting the merged bam file into 2 halves
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

while getopts ":f:h" opt; do
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

if [ ! -s reports/summary_ChIP_peaks.txt ]; then
	printf "Line\tTissue\tMark\tPeaks_in_Rep1\tPeaks_in_Rep2\tCommon_peaks\tCommon_peaks_IDR_0.05\tPeaks_in_merged\tPeaks_in_pseudo_reps\tSelected_peaks\n" > reports/summary_ChIP_peaks.txt
fi

pidsa=()
while read data line tissue mark paired ref_dir
do
	#### To merge bam files of replicates
	if [[ ${data} == "ChIP_"* ]]; then
		tmp=${data#ChIP_}
		add="_${tmp}"
	else
		add=""
	fi
	export add
	export line
	export tissue
	export mark
	export ref_dir
	export name=${line}_${tissue}_${mark}
	export input=${line}_${tissue}_Input
	export paired
	nbsample=$(ls -1 mapped/${name}_Rep*.bam | wc -l | awk '{print $1}')
	nbinput=$(ls -1 mapped/${input}_Rep*${add}.bam | wc -l | awk '{print $1}')
	if [ ${nbsample} -eq 1 ]; then
		printf "\nOnly one replicate of sample found, no merging of replicates possible\n"
		export samplerep="one"
		export inputrep="one"
	elif [ ${nbsample} -gt 2 ]; then
		printf "\nMore than two (${nbsample}) replicates found, analysis temporarily unavailable\n"
		exit 1
	elif [ -s mapped/${input}_merged${add}.bam ]; then
		printf "\nReplicates of ${input}${add} already merged\n"
		export samplerep="two"
		export inputrep="two"
	elif [ ! -s mapped/${input}_merged${add}.bam ] && [ ${nbinput} -eq 2 ]; then
		printf "\nMerging replicates of ${input}${add}\n"
		samtools merge -@ $threads mapped/temp_${input}${add}.bam mapped/${input}_Rep*${add}.bam
		samtools sort -@ $threads -o mapped/${input}_merged${add}.bam mapped/temp_${input}${add}.bam
		rm -f mapped/temp_${input}${add}.bam
		samtools index -@ $threads mapped/${input}_merged${add}.bam
		export samplerep="two"
		export inputrep="two"
	elif [ ${nbinput} -eq 1 ]; then
		printf "\nOnly one replicate of ${input}${add}\nIt will be used for all replicates\n"
		export samplerep="two"
		export inputrep="one"
	elif [ ${nbinput} -eq 0 ]; then
		printf "\nNo Input file found, cannot proceed!\n"
		exit 1
	elif [ ${nbinput} -gt 2 ]; then
		printf "\nMore than two Input files found, cannot proceed!\n"
		exit 1
	fi
	printf "\nStarting single ChIP sample analysis for $name\n"
	qsub -N ${name} -V -cwd -sync y -pe threads 2 -l m_mem_free=2G -l tmp_free=50G -j y -o logs/analysis_${name}.log <<-'EOF1' &
		#!/bin/bash
		set -e -o pipefail		
		export threads=$NSLOTS
		
		if [[ "${samplerep}" == "two" ]] && [ ! -s mapped/${name}_merged.bam ]; then
			printf "\nMerging replicates of $name\n"
			samtools merge -f -@ $threads mapped/temp_${name}.bam mapped/${name}_Rep1.bam mapped/${name}_Rep2.bam
			samtools sort -@ $threads -o mapped/${name}_merged.bam mapped/temp_${name}.bam
			rm -f mapped/temp_${name}.bam
			samtools index -@ $threads mapped/${name}_merged.bam
		fi
		if [[ "${samplerep}" == "two" ]] && [ ! -s mapped/${name}_pseudo1.bam ]; then
			printf "\nSplitting $name in two pseudo-replicates\n"
			samtools view -b -h -s 1.5 -@ $threads -U mapped/temp_${name}_pseudo2.bam -o mapped/temp_${name}_pseudo1.bam mapped/${name}_merged.bam
			samtools sort -@ $threads -o mapped/${name}_pseudo1.bam mapped/temp_${name}_pseudo1.bam
			rm -f mapped/temp_${name}_pseudo1.bam
			samtools sort -@ $threads -o mapped/${name}_pseudo2.bam mapped/temp_${name}_pseudo2.bam
			rm -f mapped/temp_${name}_pseudo2.bam
		fi
		if [[ "${samplerep}" == "one" ]] && [ ! -s mapped/${name}_pseudo1.bam ]; then
			printf "\nSplitting $name in two pseudo-replicates\n"
			samtools view -b -h -s 1.5 -@ $threads -U mapped/temp_${name}_pseudo2.bam -o mapped/temp_${name}_pseudo1.bam mapped/${name}_Rep1.bam
			samtools sort -@ $threads -o mapped/${name}_pseudo1.bam mapped/temp_${name}_pseudo1.bam
			rm -f mapped/temp_${name}_pseudo1.bam
			samtools sort -@ $threads -o mapped/${name}_pseudo2.bam mapped/temp_${name}_pseudo2.bam
			rm -f mapped/temp_${name}_pseudo2.bam
		fi
		#### To call either broad or narrow peaks if not already exisiting
		case "${mark}" in
			H3K4me1|H3K27me1|H3K27me2|H3K27me3|H3K9me1|H3K9me2|H3K9me3) export peaktype="broad";;
			H3K27ac|H3K4me3) export peaktype="narrow";;
		esac
		if [[ "${samplerep}" == "two" ]]; then
			listfiletypes=("merged" "Rep1" "Rep2" "pseudo1" "pseudo2")
		elif [[ "${samplerep}" == "one" ]]; then
			listfiletypes=("Rep1" "pseudo1" "pseudo2")
		fi
		pidsb=()
		for filetype in ${listfiletypes[@]}
		do
			export filetype
			if [[ "$inputrep" == "two" ]]; then
				case "$filetype" in
					Rep1|Rep2) 	export namefiletype=mapped/${name}_${filetype}.bam
							export inputfiletype=mapped/${input}_${filetype}${add}.bam
							export param=""
							export clean="No";;
					pseudo1|pseudo2)	export namefiletype=mapped/${name}_${filetype}.bam
									export inputfiletype=mapped/${input}_merged${add}.bam
									export param=""
									export clean="Yes";;
					merged)	export namefiletype=mapped/${name}_${filetype}.bam
						export inputfiletype=mapped/${input}_${filetype}${add}.bam
						export param="-B"
						export clean="No";;
				esac
			elif [[ "$inputrep" == "one" ]]; then
				case "$filetype" in
					Rep1) 	export namefiletype=mapped/${name}_${filetype}.bam
						export inputfiletype=mapped/${input}_${filetype}${add}.bam
						export param=""
						export clean="No";;
					Rep2)	export namefiletype=mapped/${name}_${filetype}.bam
						export inputfiletype=mapped/${input}_Rep1${add}.bam
						export param=""
						export clean="No";;
					pseudo1|pseudo2)	export namefiletype=mapped/${name}_${filetype}.bam
									export inputfiletype=mapped/${input}_Rep1${add}.bam
									export param=""
									export clean="Yes";;
					merged)	export namefiletype=mapped/${name}_${filetype}.bam
						export inputfiletype=mapped/${input}_Rep1${add}.bam
						export param="-B"
						export clean="No";;
				esac
			fi
			printf "\nStarting single ChIP sample analysis for ${name} ${filetype} using ${input}${add}\n"
			qsub -N ${name}_${filetype} -V -cwd -sync y -pe threads 10 -l m_mem_free=6G -l tmp_free=50G -j y -o logs/analysis_${name}_${filetype}.log <<-'EOF2' &
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
						printf "\nSomething is wrong with the information about peak type to call! Check usage:\n"
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
				if [[ $clean == "Yes" ]]; then
					#### To delete unnecessary bam files for pseudo-replicates
					rm -f ${namefiletype}*
				else
					#### To create bw files if not already exisiting (not for pseudo-replicates)
					if [ ! -s tracks/${name}_${filetype}.bw ]; then
						printf "\nMaking bigwig files for $namefiletype with deeptools version:\n"
						deeptools --version
						bamCompare -b1 ${namefiletype} -b2 ${inputfiletype} -o tracks/${name}_${filetype}.bw -p $threads --binSize 1 --scaleFactorsMethod "None" --normalizeUsing CPM
					else
						printf "\nBigwig file for $namefiletype already exists\n"
					fi
					#### To create fingerprint plots if not already exisiting (not for pseudo-replicates)
					if [ ! -s plots/Fingerprint_${name}_${filetype}.png ]; then
						printf "\nPlotting fingerprint for $namefiletype with deeptools version:\n"
						deeptools --version
						plotFingerprint -b ${namefiletype} ${inputfiletype} -o plots/Fingerprint_${name}_${filetype}.png -p $threads -l ${name} ${input}
					else
						printf "\nFingerprint plot for $namefiletype already exists\n"
					fi
				fi
			EOF2
			pidsb+=("$!")
		done
		printf "\nWaiting for $name single sample files to be processed\n"
		wait ${pidsb[*]}
	
		#### To get IDR analysis on biological replicates
		if [[ "${samplerep}" == "two" ]] && [ ! -s peaks/idr_${name}.${peaktype}Peak ]; then
			printf "\nDoing IDR analysis on both replicates from ${line}_${tissue}_${mark} ($peaktype peaks) with idr version:\n"
			idr --version
			idr --input-file-type ${peaktype}Peak --output-file-type ${peaktype}Peak --samples peaks/${name}_Rep1_peaks.${peaktype}Peak peaks/${name}_Rep2_peaks.${peaktype}Peak -o peaks/idr_${name}.${peaktype}Peak -l reports/idr_${name}.log --plot || true
			if [ -s peaks/idr_${name}.${peaktype}Peak.png ]; then
				mv peaks/idr_${name}.${peaktype}Peak.png plots/
			fi
		elif [ -s peaks/idr_${name}.${peaktype}Peak ]; then
			printf "\nIDR analysis already done for ${name}\n"
		fi
		#### To get the final selected peak file (peaks called in merged also present in both pseudo replicates)
		if [[ "${samplerep}" == "two" ]]; then
			awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_merged_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u > peaks/temp_${name}_merged.bed
			awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_pseudo1_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u > peaks/temp_${name}_pseudo1.bed
			awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_pseudo2_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u > peaks/temp_${name}_pseudo2.bed
			bedtools intersect -a peaks/temp_${name}_pseudo1.bed -b peaks/temp_${name}_pseudo2.bed > peaks/temp_${name}_pseudos.bed
			bedtools intersect -a peaks/temp_${name}_merged.bed -b peaks/temp_${name}_pseudos.bed -u > peaks/temp_${name}_selected.bed
			bedtools intersect -a peaks/${name}_merged_peaks.${peaktype}Peak -b peaks/temp_${name}_selected.bed -u > peaks/selected_peaks_${name}.${peaktype}Peak
		elif [[ "${samplerep}" == "one" ]]; then
			awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_Rep1_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u > peaks/temp_${name}_rep1.bed
			awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_pseudo1_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u > peaks/temp_${name}_pseudo1.bed
			awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_pseudo2_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u > peaks/temp_${name}_pseudo2.bed
			bedtools intersect -a peaks/temp_${name}_pseudo1.bed -b peaks/temp_${name}_pseudo2.bed > peaks/temp_${name}_pseudos.bed
			bedtools intersect -a peaks/temp_${name}_rep1.bed -b peaks/temp_${name}_pseudos.bed -u > peaks/temp_${name}_selected.bed
			bedtools intersect -a peaks/${name}_Rep1_peaks.${peaktype}Peak -b peaks/temp_${name}_selected.bed -u > peaks/selected_peaks_${name}.${peaktype}Peak
		fi
		printf "Getting best peak for $name\n"
		sort -k1,1 -k2,2n -k5nr peaks/selected_peaks_${name}.${peaktype}Peak | awk -v OFS="\t" '{print $1";"$2";"$3,$4,$5,$6,$7,$8,$9,$10}' | awk 'BEGIN {a=0} {b=$1; if (b!=a) print $0; a=$1}' | awk -F"[;\t]" -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | bedtools sort -g ${ref_dir}/chrom.sizes > peaks/best_peaks_${name}.bed
		
		#### To get some peaks stats for each mark
		printf "\nCalculating peak stats for ${name} in ${peaktype} peaks\n"
		if [[ "${samplerep}" == "two" ]]; then
			rep1=$(awk '{print $1,$2,$3}' peaks/${name}_Rep1_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
			rep2=$(awk '{print $1,$2,$3}' peaks/${name}_Rep2_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
			common=$(awk '{print $1,$2,$3}' peaks/idr_${name}.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
			idr=$(awk '$5>=540 {print $1,$2,$3}' peaks/idr_${name}.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
			merged=$(awk '{print $1,$2,$3}' peaks/${name}_merged_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
			pseudos=$(awk '{print $1,$2,$3}' peaks/temp_${name}_pseudos.bed | sort -k1,1 -k2,2n -u | wc -l)
			selected=$(cat peaks/temp_${name}_selected.bed | sort -k1,1 -k2,2n -u | wc -l)
			awk -v OFS="\t" -v a=$line -v b=$tissue -v c=$mark -v d=$rep1 -v e=$rep2 -v f=$common -v g=$idr -v h=$merged -v i=$pseudos -v j=$selected 'BEGIN {print a,b,c,d,e,f" ("f/d*100"%rep1;"f/e*100"%rep2)",g" ("g/f*100"%common)",h,i,j" ("j/h*100"%merged)"}' >> reports/summary_ChIP_peaks.txt
		elif [[ "${samplerep}" == "one" ]]; then
			rep1=$(awk '{print $1,$2,$3}' peaks/${name}_Rep1_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
			rep2=0
			common=0
			idr=0
			merged=$(awk '{print $1,$2,$3}' peaks/${name}_Rep1_peaks.${peaktype}Peak | sort -k1,1 -k2,2n -u | wc -l)
			pseudos=$(awk '{print $1,$2,$3}' peaks/temp_${name}_pseudos.bed | sort -k1,1 -k2,2n -u | wc -l)
			selected=$(cat peaks/temp_${name}_selected.bed | sort -k1,1 -k2,2n -u | wc -l)
			awk -v OFS="\t" -v a=$line -v b=$tissue -v c=$mark -v d=$rep1 -v e=$rep2 -v f=$common -v g=$idr -v h=$merged -v i=$pseudos -v j=$selected 'BEGIN {print a,b,c,d,e,f" ("f/d*100"%rep1;0%rep2)",g" (0%common)",h,i,j" ("j/h*100"%merged)"}' >> reports/summary_ChIP_peaks.txt
		fi
		rm -f peaks/temp_${name}*
		touch chkpts/analysis_${name}
	EOF1
	pidsa+=("$!")
done < $samplefile

printf "\nWaiting for each sample to be processed individually\n"
wait ${pidsa[*]}
printf "\nScript finished successfully!\n"

