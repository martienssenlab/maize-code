#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 20
#$ -l m_mem_free=3G
#$ -l tmp_free=20G
#$ -o TFanalysis.log
#$ -j y
#$ -N TFanalysis

usage="
##### Script for Maize code TF ChIP data analysis, used by script MaizeCode_analysis.sh for TF data
#####
##### sh MaiCode_TF_analysis.sh -f samplefile [-h]
#####	-f: samplefile containing the samples to compare and the reference directory in 6 tab-delimited columns:
##### 		Line, TF, ChIP, PE or SE, Reference directory
##### 	-h: help, returns usage
##### 
##### It merges the two replicate files, and creates pseudo-replicates by splitting the merged bam file into 2 halves
##### It calls narrow peaks with Macs2 on each biological replicate,
##### on the merged file, and on both pseudo-replicates and will keep the peaks that are present in merged and both pseudo-replicates for downstream analysis,
##### it then creates bigwig files (log2 FC vs Input) for each biological rep and the merged file, makes fingerprint plots for each,
##### makes idr analysis between biological replicates, and calculates peak stats
#####
##### Requirements: samtools, bedtools, deeptools, macs2, idr, meme, homer
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS

if [ $# -eq 0 ]; then
	printf "${usage}\n"
	exit 1
fi

while getopts ":f:h" opt; do
	case ${opt} in
		h) 	printf "${usage}\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		*)	printf "${usage}\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! ${samplefile} ]; then
	printf "Samplefile missing!\n"
	printf "${usage}\n"
	exit 1
fi

tmp1=${samplefile##*temp_}
export samplename=${tmp1%%_TF*}

if [ ! -s reports/summary_TF_peaks.txt ]; then
	printf "Line\tSample\tPeaks_in_Rep1\tPeaks_in_Rep2\tCommon_peaks\tCommon_peaks_IDR_0.05\tPeaks_in_merged\tPeaks_in_pseudo_reps\tSelected_peaks\n" > reports/summary_TF_peaks.txt
fi

pidsa=()
while read line tf chip paired ref_dir
do
	#### To merge bam files of replicates
	export line
	export tf
	export chip
	export ref_dir
	export ref=${ref_dir##*/}
	export name=${line}_${tf}
	export file=${line}_${tf}_IP
	export input=${line}_${tf}_Input
	export paired
	if [ -s mapped/${input}_merged.bam ]; then
		printf "\nReplicates of ${input} already merged\n"
		export inputrep="two"
	elif [ ! -s mapped/${input}_merged.bam ] && [ -e mapped/${input}_Rep2.bam ]; then
		printf "\nMerging replicates of ${input}\n"
		samtools merge -@ ${threads} mapped/temp_${input}.bam mapped/${input}_Rep1.bam mapped/${input}_Rep2.bam
		samtools sort -@ ${threads} -o mapped/${input}_merged.bam mapped/temp_${input}.bam
		rm -f mapped/temp_${input}.bam
		samtools index -@ ${threads} mapped/${input}_merged.bam
		export inputrep="two"
	elif [ -e mapped/${input}_Rep1.bam ] && [ ! -e mapped/${input}_Rep2.bam ]; then
		printf "\nOnly one replicate of ${input}\nIt will be used for all replicates\n"
		export inputrep="one"
	elif [ ! -e mapped/${input}_Rep1.bam ]; then
		printf "\nNo Input file found, cannot proceed!\n"
		exit 1
	fi
	printf "\nStarting single TF ChIP sample analysis for ${name}\n"
	qsub -N ${name} -V -cwd -sync y -pe threads 2 -l m_mem_free=2G -l tmp_free=50G -j y -o logs/analysis_${name}.log <<-'EOF1' &
		#!/bin/bash
		set -e -o pipefail		
		export threads=$NSLOTS
		
		if [ ! -s mapped/${name}_merged.bam ]; then
			printf "\nMerging replicates of ${name}\n"
			samtools merge -f -@ ${threads} mapped/temp_${name}.bam mapped/${file}_Rep1.bam mapped/${file}_Rep2.bam
			samtools sort -@ ${threads} -o mapped/${name}_merged.bam mapped/temp_${name}.bam
			rm -f mapped/temp_${name}.bam
			samtools index -@ ${threads} mapped/${name}_merged.bam
		fi
		if [ ! -s mapped/${name}_pseudo1.bam ]; then
			printf "\nSplitting ${name} in two pseudo-replicates\n"
			samtools view -b -h -s 1.5 -@ ${threads} -U mapped/temp_${name}_pseudo2.bam -o mapped/temp_${name}_pseudo1.bam mapped/${name}_merged.bam
			samtools sort -@ ${threads} -o mapped/${file}_pseudo1.bam mapped/temp_${name}_pseudo1.bam
			rm -f mapped/temp_${name}_pseudo1.bam
			samtools sort -@ ${threads} -o mapped/${file}_pseudo2.bam mapped/temp_${name}_pseudo2.bam
			rm -f mapped/temp_${name}_pseudo2.bam
		fi

		pidsb=()
		for filetype in merged Rep1 Rep2 pseudo1 pseudo2
		do
			export filetype
			if [[ "${inputrep}" == "two" ]]; then
				case "${filetype}" in
					Rep1|Rep2) 	export namefiletype=mapped/${file}_${filetype}.bam
							export inputfiletype=mapped/${input}_${filetype}.bam
							export param=""
							export clean="No";;
					pseudo1|pseudo2)	export namefiletype=mapped/${file}_${filetype}.bam
									export inputfiletype=mapped/${input}_merged.bam
									export param=""
									export clean="Yes";;
					merged)	export namefiletype=mapped/${name}_${filetype}.bam
						export inputfiletype=mapped/${input}_${filetype}.bam
						export param="-B"
						export clean="No";;
				esac
			elif [[ "${inputrep}" == "one" ]]; then
				case "${filetype}" in
					Rep1) 	export namefiletype=mapped/${file}_${filetype}.bam
						export inputfiletype=mapped/${input}_${filetype}.bam
						export param=""
						export clean="No";;
					Rep2)	export namefiletype=mapped/${file}_${filetype}.bam
						export inputfiletype=mapped/${input}_Rep1.bam
						export param=""
						export clean="No";;
					pseudo1|pseudo2)	export namefiletype=mapped/${file}_${filetype}.bam
									export inputfiletype=mapped/${input}_Rep1.bam
									export param=""
									export clean="Yes";;
					merged)	export namefiletype=mapped/${name}_${filetype}.bam
						export inputfiletype=mapped/${input}_Rep1.bam
						export param="-B"
						export clean="No";;
				esac
			fi
			printf "\nStarting single ChIP sample analysis for ${name} ${filetype}\n"
			qsub -N ${name}_${filetype} -V -cwd -sync y -pe threads 10 -l m_mem_free=6G -l tmp_free=50G -j y -o logs/analysis_${name}_${filetype}.log <<-'EOF2' &
				#!/bin/bash
				set -e -o pipefail
				export threads=$NSLOTS
		
				if [[ ${paired} == "PE" ]]; then
					if [ ! -s peaks/${namefiletype}_peaks.narrowPeak ]; then
						printf "\nCalling narrow peaks for PE ${namefiletype} (vs ${inputfiletype}) with macs2 version:\n"
						macs2 --version
						macs2 callpeak -t ${namefiletype} -c ${inputfiletype} -f BAMPE -g 2.2e9 ${param} -n ${name}_${filetype} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR
					elif [ -s peaks/${name}_${filetype}_peaks.narrowPeak ]; then
						printf "\nPeaks already called for ${namefiletype}\n"
					else
						printf "\nSomething is wrong with the information about peak type to call! Check usage:\n"
						printf "${usage}\n"
						exit 1
					fi
				elif [[ ${paired} == "SE" ]]; then
					if [ ! -s peaks/${name}_${filetype}_peaks.narrowPeak ]; then
						printf "\nCalling narrow peaks for SE ${namefiletype} (vs ${inputfiletype}) with macs2 version:\n"
						macs2 --version
						macs2 callpeak -t ${namefiletype} -c ${inputfiletype} -f BAM -g 2.2e9 ${param} -n ${name}_${filetype} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR --nomodel --extsize 150
					elif [ -s peaks/${name}_${filetype}_peaks.narrowPeak ]; then
						printf "\nPeaks already called for ${namefiletype}\n"
					else
						printf "\nSomething is wrong! Check usage:\n"
						printf "${usage}\n"
						exit 1
					fi
				else
					printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
					exit 1
				fi
				if [[ ${clean} == "Yes" ]]; then
					#### To delete unnecessary bam files for pseudo-replicates
					rm -f ${namefiletype}*
				else
					#### To create bw files if not already exisiting (not for pseudo-replicates)
					if [ ! -s tracks/${name}_${filetype}.bw ]; then
						printf "\nMaking bigwig files for ${namefiletype} with deeptools version:\n"
						deeptools --version
						bamCompare -b1 ${namefiletype} -b2 ${inputfiletype} -o tracks/${name}_${filetype}.bw -p $threads --binSize 1 --scaleFactorsMethod "None" --normalizeUsing CPM
					else
						printf "\nBigwig file for ${namefiletype} already exists\n"
					fi
					#### To create fingerprint plots if not already exisiting (not for pseudo-replicates)
					if [ ! -s plots/Fingerprint_${name}_${filetype}.png ]; then
						printf "\nPlotting fingerprint for ${namefiletype} with deeptools version:\n"
						deeptools --version
						plotFingerprint -b ${namefiletype} ${inputfiletype} -o plots/Fingerprint_${name}_${filetype}.png -p ${threads} -l ${name} ${input}
					else
						printf "\nFingerprint plot for ${namefiletype} already exists\n"
					fi
				fi
			EOF2
			pidsb+=("$!")
		done
		printf "\nWaiting for ${name} single sample files to be processed\n"
		wait ${pidsb[*]}
	
		#### To get IDR analysis on biological replicates
		if [ ! -s peaks/idr_${name}.narrowPeak ]; then
			printf "\nDoing IDR analysis on both replicates from ${line}_${tf}_${chip} with idr version:\n"
			idr --version
			idr --input-file-type narrowPeak --output-file-type narrowPeak --samples peaks/${name}_Rep1_peaks.narrowPeak peaks/${name}_Rep2_peaks.narrowPeak -o peaks/idr_${name}.narrowPeak -l reports/idr_${name}.log --plot || true
			if [ -s peaks/idr_${name}.narrowPeak.png ]; then
				mv peaks/idr_${name}.narrowPeak.png plots/
			fi
		else
			printf "\nIDR analysis already done for ${name}\n"
		fi
		#### To get the final selected peak file (peaks called in merged also present in both pseudo replicates)
		awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_merged_peaks.narrowPeak | sort -k1,1 -k2,2n -u > peaks/temp_${name}_merged.bed
		awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_pseudo1_peaks.narrowPeak | sort -k1,1 -k2,2n -u > peaks/temp_${name}_pseudo1.bed
		awk -v OFS="\t" '{print $1,$2,$3}' peaks/${name}_pseudo2_peaks.narrowPeak | sort -k1,1 -k2,2n -u > peaks/temp_${name}_pseudo2.bed
		bedtools intersect -a peaks/temp_${name}_pseudo1.bed -b peaks/temp_${name}_pseudo2.bed > peaks/temp_${name}_pseudos.bed
		bedtools intersect -a peaks/temp_${name}_merged.bed -b peaks/temp_${name}_pseudos.bed -u > peaks/temp_${name}_selected.bed
		bedtools intersect -a peaks/${name}_merged_peaks.narrowPeak -b peaks/temp_${name}_selected.bed -u > peaks/selected_peaks_${name}.narrowPeak
		printf "Getting best peak for ${name}\n"
		sort -k1,1 -k2,2n -k5nr peaks/selected_peaks_${name}.narrowPeak | awk -v OFS="\t" '{print $1";"$2";"$3,$4,$5,$6,$7,$8,$9,$10}' | awk 'BEGIN {a=0} {b=$1; if (b!=a) print $0; a=$1}' | awk -F"[;\t]" -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > peaks/best_peaks_${name}.bed

		#### To get some peaks stats for each TF
		printf "\nCalculating peak stats for ${name} in narrow peaks\n"
		rep1=$(awk '{print $1,$2,$3}' peaks/${name}_Rep1_peaks.narrowPeak | sort -k1,1 -k2,2n -u | wc -l)
		rep2=$(awk '{print $1,$2,$3}' peaks/${name}_Rep2_peaks.narrowPeak | sort -k1,1 -k2,2n -u | wc -l)
		common=$(awk '{print $1,$2,$3}' peaks/idr_${name}.narrowPeak | sort -k1,1 -k2,2n -u | wc -l)
		idr=$(awk '$5>=540 {print $1,$2,$3}' peaks/idr_${name}.narrowPeak | sort -k1,1 -k2,2n -u | wc -l)
		merged=$(awk '{print $1,$2,$3}' peaks/${name}_merged_peaks.narrowPeak | sort -k1,1 -k2,2n -u | wc -l)
		pseudos=$(awk '{print $1,$2,$3}' peaks/temp_${name}_pseudos.bed | sort -k1,1 -k2,2n -u | wc -l)
		selected=$(cat peaks/temp_${name}_selected.bed | sort -k1,1 -k2,2n -u | wc -l)
		if [[ ${rep1} -gt 0 ]] && [[ ${rep2} -gt 0 ]] && [[ ${common} -gt 0 ]] && [[ ${merged} -gt 0 ]]; then
			awk -v OFS="\t" -v a=${line} -v b=${tf} -v d=${rep1} -v e=${rep2} -v f=${common} -v g=${idr} -v h=${merged} -v i=${pseudos} -v j=${selected} 'BEGIN {print a,b,d,e,f" ("f/d*100"%rep1;"f/e*100"%rep2)",g" ("g/f*100"%common)",h,i,j" ("j/h*100"%merged)"}' >> reports/summary_TF_peaks.txt
		else
			awk -v OFS="\t" -v a=${line} -v b=${tf} -v d=${rep1} -v e=${rep2} -v f=${common} -v g=${idr} -v h=${merged} -v i=${pseudos} -v j=${selected} 'BEGIN {print a,b,d,e,f,g,h,i,j}' >> reports/summary_TF_peaks.txt		
		fi
		rm -f peaks/temp_${name}*
		
		#### To find motifs in different peak sets:
		if [ -s ${ref_dir}/temp_${ref}.fa ]; then
			fasta=${ref_dir}/temp_${ref}.fa
		elif [ -s ${ref_dir}/*.fa.gz ]; then
			fa_file=$(ls ${ref_dir}/*.fa.gz)
			pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${ref}.fa
			fasta=${ref_dir}/temp_${ref}.fa
		elif [ -s ${ref_dir}/*.fa ]; then
			fa_file=$(ls ${ref_dir}/*.fa)
			fasta=${fa_file}
		elif [ -s ${ref_dir}/*.fasta.gz ]; then
			fa_file=$(ls ${ref_dir}/*.fasta.gz)
			pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${ref}.fa
			fasta=${ref_dir}/temp_${ref}.fa
		elif [ -s ${ref_dir}/*.fasta ]; then
			fa_file=$(ls ${ref_dir}/*.fasta)
			fasta=${fa_file}
		fi

		if [ ! -d motifs/${name}/meme ] && [ -s peaks/best_peaks_${name}.bed ]; then
			#### v1="selected" peaks (best peaks from selected, i.e. in merged and both pseudo reps) with MEME		
			printf "\nGetting peak fasta sequences for ${name} meme v1\n"
			awk -v OFS="\t" '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {a=$2+$10; print $1,a-50,a+50,$4}' peaks/best_peaks_${name}.bed > peaks/selected_motifs_${name}.bed
			bedtools getfasta -name -fi ${fasta} -bed peaks/selected_motifs_${name}.bed > peaks/selected_sequences_${name}.fa
			printf "\nGetting motifs for ${name} with meme\n"
			meme-chip -oc motifs/${name}/meme -meme-p ${threads} -meme-nmotifs 10 -streme-nmotifs 10 peaks/selected_sequences_${name}.fa
			if [ -s motifs/${name}/meme/combined.meme ]; then
				printf "\nLooking for similar motifs in JASPAR database with tomtom\n"
				tomtom -oc motifs/${name}/tomtom motifs/${name}/meme/combined.meme motifs/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt
			fi
		fi
		if [ ! -d motifs/${name}/meme2 ] && [ -s peaks/idr_${name}.narrowPeak ]; then
			#### v2="replicated" peaks (peaks in both biological reps, i.e all peaks in idr) with MEME
			printf "\nGetting peak fasta sequences for ${name} meme v2\n"
			awk -v OFS="\t" '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {a=$2+$10; print $1,a-50,a+50}' peaks/idr_${name}.narrowPeak > peaks/selected_motifs_${name}.bed
			bedtools getfasta -name -fi ${fasta} -bed peaks/selected_motifs_${name}.bed > peaks/selected_sequences_${name}.fa
			printf "\nGetting motifs for ${name} with meme\n"
			meme-chip -oc motifs/${name}/meme2 -meme-p ${threads} -meme-nmotifs 10 -streme-nmotifs 10 peaks/selected_sequences_${name}.fa
			if [ -s motifs/${name}/meme2/combined.meme ]; then
				printf "\nLooking for similar motifs in JASPAR database with tomtom\n"
				tomtom -oc motifs/${name}/tomtom2 motifs/${name}/meme2/combined.meme motifs/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt
			fi
		fi
		
		#### Needs a file of masked regions. Can be found online, or created with Repeat masker. 
		#### We could consider having a helper script, or documentation on how to create it.
		#### For now just copying the file generated from repeat masker with the following commands:
		# RepeatMasker -pa $((threads/4)) -lib MTEC/maizeTE02052020 -gccalc -xsmall -q -no_is -norna -div 40 -nolow -html -gff -poly -dir B73_v4 originals/B73_v4.fa >> logs/B73_v4_repeatmasker_01.log 2>&1
		# awk -v OFS="\t" '$1 !~ /^#/ {print $1,$4,$5,$9,".","."}' ~/nlsas/Genomes/RepeatMasker/B73_v4/B73_v4.fa.out.gff > tracks/B73_v4_masked_regions.bed
		
		#### For homer, it requires to load the masked genome first with the following command:
		# loadGenome.pl -name B73_v4 -org null -fasta ~/nlsas/Genomes/RepeatMasker/B73_v4/B73_v4.fa.masked -gtf ~/nlsas/Genomes/Zea_mays/B73_v4/temp_B73_v4.gtf -promoters "B73_v4-p"
		
		if [[ ${ref} == "B73_v4" ]]; then
			if [ ! -d motifs/${name}/meme3 ] && [ -s peaks/best_peaks_${name}.bed ]; then
				#### v3="selected" peaks (best peaks from selected, i.e. in merged and both pseudo reps) without masked sequences with MEME
				printf "\nGetting peak fasta sequences for ${name} meme v3\n"
				awk -v OFS="\t" '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {a=$2+$10; print $1,a-50,a+50,$4}' peaks/best_peaks_${name}.bed > peaks/selected_motifs_${name}.bed
				bedtools intersect -v -wa -a peaks/selected_motifs_${name}.bed -b tracks/${ref}_masked_regions.bed > peaks/masked_selected_motifs_${name}.bed
				bedtools getfasta -name -fi ${fasta} -bed peaks/masked_selected_motifs_${name}.bed > peaks/selected_sequences_${name}.fa
				printf "\nGetting motifs for ${name} with meme\n"
				meme-chip -oc motifs/${name}/meme3 -meme-p ${threads} -meme-nmotifs 10 -streme-nmotifs 10 peaks/selected_sequences_${name}.fa
				if [ -s motifs/${name}/meme3/combined.meme ]; then
					printf "\nLooking for similar motifs in JASPAR database with tomtom\n"
					tomtom -oc motifs/${name}/tomtom3 motifs/${name}/meme3/combined.meme motifs/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt
				fi
			fi
			if [ ! -d motifs/${name}/meme4 ] && [ -s peaks/idr_${name}.narrowPeak ]; then
				#### v4="replicated" peaks (peaks in both biological reps, i.e all peaks in idr) without masked sequences with MEME
				printf "\nGetting peak fasta sequences for ${name} meme v4\n"
				awk -v OFS="\t" '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {a=$2+$10; print $1,a-50,a+50}' peaks/idr_${name}.narrowPeak > peaks/selected_motifs_${name}.bed
				bedtools intersect -v -wa -a peaks/selected_motifs_${name}.bed -b tracks/${ref}_masked_regions.bed > peaks/masked_selected_motifs_${name}.bed
				bedtools getfasta -name -fi ${fasta} -bed peaks/masked_selected_motifs_${name}.bed > peaks/selected_sequences_${name}.fa
				printf "\nGetting motifs for ${name} with meme\n"
				meme-chip -oc motifs/${name}/meme4 -meme-p ${threads} -meme-nmotifs 10 -streme-nmotifs 10 peaks/selected_sequences_${name}.fa
				if [ -s motifs/${name}/meme4/combined.meme ]; then
					printf "\nLooking for similar motifs in JASPAR database with tomtom\n"
					tomtom -oc motifs/${name}/tomtom4 motifs/${name}/meme4/combined.meme motifs/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt
				fi
			fi
			if [ ! -d motifs/${name}/homer ] && [ -s peaks/best_peaks_${name}.bed ]; then
				#### v5="selected" peaks (best peaks from selected, i.e. in merged and both pseudo reps) with HOMER		
				printf "\nGetting motifs for ${name} with HOMER\n"
				findMotifsGenome.pl peaks/best_peaks_${name}.bed B73_v4 motifs/${name}/homer -len 6,8,10 -size given -p ${thread}
			fi
		fi
		
##### The problem with this is that the best motif(s) is not necessarily the first one(s), so an alternative approach should be developped
##### An alternative is to have people look at the results and choose their favorite motifs and extract the corresponding peaks with a helper script, and then repeat the analysis.
#		if [ -e motifs/peaks_with_motifs_${name}_meme2.txt ]; then
#			rm -f motifs/peaks_with_motifs_${name}_meme2.txt
#		fi
#		for num in 1 2 3 ### numer of the motif(s) to select
#		do
#			awk -v OFS="\t" -v n=$num 'NR > 1 && $1 !~ /^#/ {print $3,$4,$5,"motif_"n}' motifs/${name}/meme2/fimo_out_${num}/fimo.tsv | sort -k1,1 -k2,2n > motifs/temp_motifs_${name}.bed
#			if [ -s motifs/temp_motifs_${name}.bed ]; then
#				printf "Extracting peaks containing motif number ${num}\n"
#				sort -k1,1 -k2,2n peaks/idr_${name}.narrowPeak > peaks/temp_sorted_idr_${name}.bed
#				bedtools intersect -loj -a peaks/temp_sorted_idr_${name}.bed -b motifs/temp_motifs_${name}.bed >> motifs/peaks_with_motifs_${name}_meme2.txt
#				rm -f motifs/temp_motifs_${name}.bed peaks/temp_sorted_idr_${name}.bed
#			fi
#		done
		
		touch chkpts/analysis_${name}
	EOF1
	pidsa+=("$!")
done < $samplefile

printf "\nWaiting for each sample to be processed individually\n"
wait ${pidsa[*]}

if [ -e ${ref_dir}/temp_${ref}.fa ]; then
	rm -f ${ref_dir}/temp_${ref}.fa
fi
printf "\nScript finished successfully!\n"

