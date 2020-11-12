#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 1
#$ -l m_mem_free=1G
#$ -l tmp_free=1G
#$ -o random.log
#$ -j y
#$ -N random

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS
printf "Original threads: $threads\n"

for name in B73_endosperm_RNAseq 
do
	#### To get some stats for each RNAseq sample
	printf "\nCalculating gene expession stats for ${name}\n"
	qsub -N test -V -cwd -sync y -pe threads 2 -l m_mem_free=2G -l tmp_free=5G -j y -o random.log <<-EOF
		#!/bin/bash
		set -e -o pipefail
		
		export threads=\$NSLOTS
		printf "\nCalculating gene expession stats for ${name}\n"
		totrep1=$(awk '$1 !~ /^N/ {n+=$2} END {print n}' mapped/map_${name}_Rep1_ReadsPerGene.out.tab)
		printf "totrep1: \${totrep1}\n\n"
		awk -v OFS="\t" -v t=\${totrep1} '$1 !~ /^N/ {n=$2/t; print $1,n}' mapped/map_${name}_Rep1_ReadsPerGene.out.tab > mapped/temp_genes_${name}_rep1.txt
		rep1=$(awk -v OFS="\t" '{if ($2=0) a+=1; else if ($2<1) b+=1; else c+=1} END {printf "a\tb\tc"}' mapped/temp_genes_${name}_rep1.txt)
		totrep2=$(awk '$1 !~ /^N_/ {n+=$2} END {print n}' mapped/map_${name}_Rep2_ReadsPerGene.out.tab)
		awk -v OFS="\t" -v t=\$totrep2 '$1 !~ /^N_/ {n=$2/t; print $1,n}' mapped/map_${name}_Rep2_ReadsPerGene.out.tab > mapped/temp_genes_${name}_rep2.txt
		rep2=$(awk -v OFS="\t" '{if ($2=0) a+=1; else if ($2<1) b+=1; else c+=1} END {printf "a\tb\tc"}' mapped/temp_genes_${name}_rep2.txt)
		while read gene value
		do
			awk -v OFS="\t" -v g=\$gene -v r=\$value '$1==g {m=($2+r)/2; print $1,m}' mapped/temp_genes_${name}_rep1.txt > mapped/temp_genes_${name}_mean.txt
		done < mapped/temp_genes_${name}_rep1.txt
		mean=$(awk -v OFS="\t" '{if ($2=0) a+=1; else if ($2<1) b+=1; else c+=1} END {printf "a\tb\tc"}' mapped/temp_genes_${name}_mean.txt)
		awk -v OFS="\t" -v l=$line -v t=$tissue -v r=$rnatype -v n=${annotated_gene_number} -v a=\$rep1 -v b=\$rep2 -v c=\$mean 'BEGIN {printf l,t,r,n,a,b,c}' >> TSS/summary_expression_${samplename}.txt
		rm -f mapped/temp_genes_${name}_*
EOF
done
