#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=2G
#$ -l tmp_free=10G
#$ -o RNAanalysis.log
#$ -j y
#$ -N RNAanalysis

usage="
##### Script for Maize code RNA data analysis, used by script MaizeCode_analysis.sh for RNA data
#####
##### sh MaiCode_RNA_analysis.sh -f samplefile [-h]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Mark, PE or SE, Reference genome directory
##### 	-h: help, returns usage
##### 
##### It merges the two replicate files and creates stranded bigwig files for each sample
##### For RNAseq samples, it creates a list of expressed genes in each biological replicates and select expressed genes (average cpm in biological replicates >1)
##### For RAMPAGE, it calls peaks to identify TSS in each biological replicate and makes IDR analysis between replicate and makes stats
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
export samplename=${tmp1%%_RNA*}

if [ ! -s reports/summary_RAMPAGE_tss.txt ]; then
	printf "Line\tTissue\tType\tTotal_annotated_genes\tTSS_in_rep1\tTSS_in_Rep2\tCommon_TSS\tCommon_TSS_IDR<=0.05\n" > reports/summary_RAMPAGE_tss.txt
fi
if [ ! -s reports/summary_gene_expression.txt ]; then
	printf "Line\tTissue\tType\tTotal_annotated_genes\tNot_expressed_in_Rep1\tLow_expression_in_Rep1(<1cpm)\tHigh_expression_in_Rep1(>1cpm)\tNot_expressed_in_Rep2\tLow_expression_in_Rep2(<1cpm)\tHigh_expression_in_Rep2(>1cpm)\tNo_mean\tLow_mean(<1cpm)\tHigh_mean(>1cpm)\n" > reports/summary_gene_expression.txt
fi

pids=()
while read line tissue rnatype paired ref_dir
do
	export line
	export tissue
	export rnatype
	export name=${line}_${tissue}_${rnatype}
	case "$rnatype" in
		RNAseq)	export param_bg="--outWigType bedGraph"
				export strandedness="reverse";;
		RAMPAGE)	export param_bg="--outWigType bedGraph read1_5p"
				export strandedness="forward";;
	esac
	printf "\nStarting single RNA sample analysis for $name\n"	
	export ref_dir=${ref_dir}
	export ref=${ref_dir##*/}
	export annotated_gene_number=$(cat tracks/${ref}_all_genes.bed | wc -l)
	qsub -N ${name} -V -cwd -sync y -pe threads 20 -l m_mem_free=5G -l tmp_free=50G -j y -o logs/analysis_${name}.log <<-'EOF' &
		#!/bin/bash
		set -e -o pipefail
		
		printf "\n\n"
		date
		printf "\n"
		
		export threads=$NSLOTS
		
		#### To merge bam files of replicates
		if [ ! -s mapped/${name}_merged.bam ]; then
			printf "\nMerging replicates of $name\n"
			if [ -e mapped/mrkdup_${name}_Rep1_Processed.out.bam ]; then
				samtools merge -@ $threads mapped/temp_${name}.bam mapped/mrkdup_${name}_Rep*_Processed.out.bam
			else
				samtools merge -@ $threads mapped/temp_${name}.bam mapped/map_${name}_Rep*_Aligned.sortedByCoord.out.bam
			fi
			samtools sort -@ $threads -o mapped/${name}_merged.bam mapped/temp_${name}.bam
			rm -f mapped/temp_${name}.bam
			samtools index -@ $threads mapped/${name}_merged.bam
		fi
		#### To make bw files of merged samples if not already existing
		if [ ! -s tracks/${name}_merged_minus.bw ]; then
			### Making BedGraph files
			printf "\nMaking bedGraph files\n"
			STAR --runMode inputAlignmentsFromBAM --inputBAMfile mapped/${name}_merged.bam --outWigStrand Stranded ${param_bg} --outFileNamePrefix tracks/bg_${name}_merged_
			### Converting to bigwig files
			printf "\nConverting bedGraphs to bigWigs\n"
			bedSort tracks/bg_${name}_merged_Signal.UniqueMultiple.str1.out.bg tracks/${name}_merged_Signal.sorted.UniqueMultiple.str1.out.bg
			bedSort tracks/bg_${name}_merged_Signal.Unique.str1.out.bg tracks/${name}_merged_Signal.sorted.Unique.str1.out.bg
			bedSort tracks/bg_${name}_merged_Signal.UniqueMultiple.str2.out.bg tracks/${name}_merged_Signal.sorted.UniqueMultiple.str2.out.bg
			bedSort tracks/bg_${name}_merged_Signal.Unique.str2.out.bg tracks/${name}_merged_Signal.sorted.Unique.str2.out.bg
			if [[ $strandedness == "forward" ]]; then
				bedGraphToBigWig tracks/${name}_merged_Signal.sorted.UniqueMultiple.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_merged_plus.bw
				bedGraphToBigWig tracks/${name}_merged_Signal.sorted.Unique.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_merged_unique_plus.bw
				bedGraphToBigWig tracks/${name}_merged_Signal.sorted.UniqueMultiple.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_merged_minus.bw
				bedGraphToBigWig tracks/${name}_merged_Signal.sorted.Unique.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_merged_unique_minus.bw
			elif [[ $strandedness == "reverse" ]]; then
				bedGraphToBigWig tracks/${name}_merged_Signal.sorted.UniqueMultiple.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_merged_minus.bw
				bedGraphToBigWig tracks/${name}_merged_Signal.sorted.Unique.str1.out.bg ${ref_dir}/chrom.sizes tracks/${name}_merged_unique_minus.bw
				bedGraphToBigWig tracks/${name}_merged_Signal.sorted.UniqueMultiple.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_merged_plus.bw
				bedGraphToBigWig tracks/${name}_merged_Signal.sorted.Unique.str2.out.bg ${ref_dir}/chrom.sizes tracks/${name}_merged_unique_plus.bw
			else
				printf "\nStrandedness of data unknown! Tracks could not be created and script failed!\n"
				exit 1
			fi
			### Moving Log files to report folder
			mv tracks/*${name}_merged*Log* reports/
			### Cleaning up
			rm -f tracks/*${name}_merged_Signal*
		else
			printf "\nBigwig files for $name already exists\n"
		fi
		
		#### Performing peak calling on RAMPAGE data to identify TSS and check biological variation with IDR
		if [[ $rnatype == "RAMPAGE" ]]; then
			#### To call TSS/peaks on each biological replicate
			for rep in Rep1 Rep2
			do
				if [ ! -s TSS/${name}_${rep}_peaks.narrowPeak ] && [ -s mapped/mrkdup_${line}_${tissue}_RNAseq_${rep}_Processed.out.bam ]; then
					namefile=mapped/mrkdup_${name}_${rep}_Processed.out.bam
					controlfile=mapped/mrkdup_${line}_${tissue}_RNAseq_${rep}_Processed.out.bam
					printf "\nCalling peaks for TSS on ${namefile} vs ${controlfile} with macs2 version:\n"
					macs2 --version
					macs2 callpeak -t ${namefile} -c ${controlfile} -f BAM -g 2.2e9 -B -n ${name}_${rep} --keep-dup "all" --call-summits --outdir TSS/ --tempdir $TMPDIR --nomodel --extsize 100
				else
					printf "\nPeak already called for $name $rep or no corresponding RNAseq control file found\n"
				fi
			done
			#### To get IDR analysis on biological replicates
			if [ ! -s TSS/idr_${name}.narrowPeak ] && [ -s TSS/${name}_Rep1_peaks.narrowPeak ] && [ -s TSS/${name}_Rep2_peaks.narrowPeak ]; then
				printf "\nGoing IDR analysis on both replicates from ${name} with idr version:\n"
				idr --version
				idr --input-file-type narrowPeak --samples TSS/${name}_Rep1_peaks.narrowPeak TSS/${name}_Rep2_peaks.narrowPeak -o TSS/idr_${name}.narrowPeak -l reports/idr_${name}.log --plot || true
				mv TSS/idr_${name}.narrowPeak.png plots/idr_${name}.png
				
				#### To get some tss stats for each RAMPAGE sample
				printf "\nCalculating peak/TSS stats for ${name}\n"
				rep1=$(awk '{print $1,$2,$3}' TSS/${name}_Rep1_peaks.narrowPeak | sort -k1,1 -k2,2n -u | wc -l)
				rep2=$(awk '{print $1,$2,$3}' TSS/${name}_Rep2_peaks.narrowPeak | sort -k1,1 -k2,2n -u | wc -l)
				common=$(awk '{print $1,$2,$3}' TSS/idr_${name}.narrowPeak | sort -k1,1 -k2,2n -u | wc -l)
				idr=$(awk '$5>=540 {print $1,$2,$3}' TSS/idr_${name}.narrowPeak | sort -k1,1 -k2,2n -u | wc -l)
				awk -v OFS="\t" -v a=$line -v b=$tissue -v c=$rnatype -v n=${annotated_gene_number} -v d=$rep1 -v e=$rep2 -v f=$common -v g=$idr 'BEGIN {print a,b,c,n,d,e,f" ("f/d*100"%rep1;"f/e*100"%rep2)",g" ("g/f*100"%common)"}' >> reports/summary_RAMPAGE_tss.txt
			else
				printf "\nIDR analysis already done for ${name} or no peak files found\n"
			fi
		#### To get some stats on gene expression levels for RNAseq samples
		elif [[ $rnatype == "RNAseq" ]]; then
			printf "\nCalculating gene expession stats for ${name}\n"
			totrep1=$(awk '$1 !~ /^N/ {n+=$2} END {print n}' mapped/map_${name}_Rep1_ReadsPerGene.out.tab)
			awk -v OFS="\t" -v t=${totrep1} '$1 !~ /^N/ {n=$2/t*1000000; print $1,n}' mapped/map_${name}_Rep1_ReadsPerGene.out.tab > mapped/temp_genes_${name}_rep1.txt
			totrep2=$(awk '$1 !~ /^N_/ {n+=$2} END {print n}' mapped/map_${name}_Rep2_ReadsPerGene.out.tab)
			awk -v OFS="\t" -v t=$totrep2 '$1 !~ /^N_/ {n=$2/t*1000000; print $1,n}' mapped/map_${name}_Rep2_ReadsPerGene.out.tab > mapped/temp_genes_${name}_rep2.txt
			paste mapped/temp_genes_${name}_rep1.txt mapped/temp_genes_${name}_rep2.txt | awk -v OFS="\t" '{m=($2+$4)/2; print $1,$2,$4,m}' > mapped/temp_genes_${name}_mean.txt
			awk -v OFS="\t" -v l=$line -v t=$tissue -v r=$rnatype -v n=${annotated_gene_number} '{if ($2==0) a+=1; else if ($2<=1) b+=1; else c+=1; if ($3==0) d+=1; else if ($3<=1) e+=1; else f+=1; if ($4==0) g+=1; else if ($4<=1) h+=1; else i+=1} END {print l,t,r,n,a,b,c,d,e,f,g,h,i}' mapped/temp_genes_${name}_mean.txt >> reports/summary_gene_expression.txt
			rm -f mapped/temp_genes_${name}*
		fi
		printf "\nAnalysis finished for $name\n"
		touch chkpts/analysis_${name}
	EOF
	pids+=("$!")
done < $samplefile

printf "\nWaiting for samples to be processed individually\n"
wait ${pids[*]}

printf "\nScript finished successfully!\n"
