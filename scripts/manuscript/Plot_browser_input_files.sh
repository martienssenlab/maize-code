#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 12
#$ -l m_mem_free=2G
#$ -l tmp_free=2G
#$ -o plotbrowser.log
#$ -j y
#$ -N browser

set -e -o pipefail

export threads=$NSLOTS

printf "\n\n"
date
printf "\n"

################################################################################
### To make browser shots of selected regions (here top 20 distal peaks in ears)
### run: qsub Plot_browser_input_files.sh -l line -r regions -f samplefile [-s]
### Samplefile should have the following columns in tab separated file: samplename group pathtobw backcolor trackcolor fillcolor1 fillcolor2
### where: 	- samplename is the name of the sample
###			- group is a unique identifier (number or text) to scale all samples with the same identifier together
###				samplename_group should be unique!
###			- pathtobw is the whole path and name of the bigwig file
###			- backcolor is the color of the background for the sample name
### 		- trackcolor is the color of the line for the sample
### 		- fillcolor1 is the color of negative values for the sample
### 		- fillcolor2 is the color of positive values for the sample
###			- invert is either "no" or "invert" to tell whether to invert the bigwig file (for stranded data)
### regions is a file of loci with minimum two columns chr:start:end and ID, and potentially highligh regions (starts, widths), 
### -s can set if the regions should be fixed or extended to all genes intersecting the loci (default=fixed)

if [ ! -d data ]; then
	mkdir data
fi

if [ ! -d plots ]; then
	mkdir plots
fi

while getopts "l:r:f:s" opt; do
	case $opt in
		l) 	export line=${OPTARG};;
		r)	export loci=${OPTARG};;
		f)	export samplefile=${OPTARG};;
		s)	export extend="Yes";;
		*)	printf "\nWrong argument: must be -l line -r loci -f samplefile\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

case "$line" in
	W22)	ref="W22_v2"
			gff="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/W22_v2_no_chr.gff"
			path="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode";;
	B73)	ref="B73_v5"
			gff="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/B73_v5_no_chr.gff"
			path="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode";;
	TIL11)	ref="TIL11_chrs"
			gff="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/TIL11_no_chr.gff"
			path="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode";;
	NC350)	ref="NC350_NAM"
			gff="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/NC350_NAM_no_chr.gff"
			path="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode";;
	B73_v4)	ref="B73_v4"
			gff="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/B73_v4_no_chr.gff"
			path="/grid/martienssen/home/jcahn/nlsas/projects/maizecode_v4";;
	Colcen)	ref="ColCEN"
			gff="/grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis/Col_CEN/ColCEN_GENES_Araport11.gff3"
			path="/grid/martienssen/data_nlsas_norepl/jcahn/projects/pipeline";;
esac
ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/${ref}"

awk '$1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/' ${path}/combined/TSS/${ref}_all_tes.bed > data/${ref}_all_tes.bed
if [[ ${ref} == "B73_v4" ]]; then
	awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,".",$7,$9}' ${gff} | awk -v OFS="\t" -F"[:;]" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$7,$4,$5}' > data/${ref}_all_genes.bed
else
	awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,".",$7,$9}' ${gff} | awk -v OFS="\t" -F"[=;]" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$7,$4,$5}' > data/${ref}_all_genes.bed
fi

filelist=()
while read sample group bw backcolor trackcolor fillcolor1 fillcolor2 invert
do
	name="${sample}_${group}"
	filelist+=("${bw}")
	labellist+=("${name}")
	printf "${name} added to filelist\n"
done < ${samplefile}
	
while read locus ID bins htstart htwidth
do
	printf "Preparing files for browser on ${ID}\n"
	rm -f data/genes_in_locus.*
	printf "${locus}\n" | awk -F"[:]" -v OFS="\t" '{print $1,$2,$3}' > data/locus_${ID}.bed
	chr=$(awk '{print $1}' data/locus_${ID}.bed)
	start=$(awk '{print $2}' data/locus_${ID}.bed)
	end=$(awk '{print $3}' data/locus_${ID}.bed)
	bedtools intersect -a data/${ref}_all_genes.bed -b data/locus_${ID}.bed | awk '{print $4}' > data/genes_in_locus_${ID}.txt
	if [ -s data/genes_in_locus_${ID}.txt ] && [[ "${extend}" == "Yes" ]]; then
		grep -f data/genes_in_locus_${ID}.txt ${gff} | awk -v OFS="\t" '{if ($7!="+" && $7!="-") $7="*"; print $0}' > data/genes_in_locus_${ID}.gff
		region=$(awk -v OFS=":" -v s=${start} -v e=${end} '{if (NR==1) {c=$1; a=$4-1;}} END {b=$5; if (a<s) m=a; else m=s; if (b>e) n=b; else n=e; print c,m,n}' data/genes_in_locus_${ID}.gff)
	elif [ -s data/genes_in_locus_${ID}.txt ]; then
		bedtools intersect -wa -a ${gff} -b data/locus_${ID}.bed | awk -v OFS="\t" '{if ($7!="+" && $7!="-") $7="*"; print $0}' > data/genes_in_locus_${ID}.gff
		region="${locus}"
	else
		region="${locus}"
	fi
	bedtools intersect -a data/${ref}_all_tes.bed -b data/locus_${ID}.bed | awk -v OFS="\t" '{if ($6!="+" && $6!="-") $6="*"; print $0}' > data/tes_in_locus_${ID}.bed
	printf "Testing if all bw have data on the chromosome\n"
	filelist2=()
	i=0
	for file in ${filelist[@]}
	do
		bigWigToBedGraph -chrom=${chr} -start=${start} -end=${end} ${file} data/${labellist[i]}_empty.bg
		if [ -s data/${labellist[i]}_empty.bg ]; then
			printf "${labellist[i]} has data on ${chr}\n"
			filelist2+=("${file}")
		else
			printf "${labellist[i]} is empty on ${chr}\n"
			grep "${chr}" ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{print $1,"1",$2,"0"}' > data/${labellist[i]}_empty.bg
			bedGraphToBigWig data/${labellist[i]}_empty.bg ${ref_dir}/chrom.sizes data/${labellist[i]}_empty.bw
			filelist2+=("data/${labellist[i]}_empty.bw")
		fi
		i=$((i+1))
		rm -f data/${labellist[i]}_empty.bg
	done
	printf "Summarize bigwigs in binsize of ${bins} bp on ${region}\n"
	multiBigwigSummary bins -b ${filelist2[@]} -l ${labellist[@]} -r ${region} -p ${threads} -bs=${bins} -out data/${ID}.npz --outRawCounts data/locus_${ID}.tab
	rm -f data/${ID}.npz

	while read sample group bw backcolor trackcolor fillcolor1 fillcolor2 invert
	do
		name="${sample}_${group}"
		printf "Making bw for ${name}\n"
		col=($(awk -v ORS=" " -v t=${name} 'NR==1 {for(i=1;i<=NF;i++) if ($i~t) print i}' data/locus_${ID}.tab))
		if [[ ${invert} == "invert" ]]; then
			awk -v OFS="\t" -v a=${col} 'NR>1 {if ($a == "nan") b=0; else b=-$a; print $1,$2,$3,b}' data/locus_${ID}.tab | bedtools sort -g ${ref_dir}/chrom.sizes > data/${ID}_${name}.bedGraph
		else
			awk -v OFS="\t" -v a=${col} 'NR>1 {if ($a == "nan") b=0; else b=$a; print $1,$2,$3,b}' data/locus_${ID}.tab | bedtools sort -g ${ref_dir}/chrom.sizes > data/${ID}_${name}.bedGraph
		fi
		bedGraphToBigWig data/${ID}_${name}.bedGraph ${ref_dir}/chrom.sizes data/${name}_locus.bw
	done < ${samplefile}
	
	printf "Name\tGroup\tBackcolor\tTrackcolor\tFillcolor1\tFillcolor2\tYmin\tYmax\n" > data/filenames_${ID}.txt
	while read sample group bw backcolor trackcolor fillcolor1 fillcolor2 invert lims
	do
		if [[ ${lims} != "" ]]; then
			ylimmin=${lims%%;*}
			ylimmax=${lims##*;}
		else
			ylimmin=$(cat data/${ID}_*${group}.bedGraph | awk 'BEGIN {a=9999} {if ($4<a) a=$4;} END {if (a<0) b=a*1.2; else b=a*0.8; print b}' )
			ylimmax=$(cat data/${ID}_*${group}.bedGraph | awk 'BEGIN {a=-9999} {if ($4>a) a=$4;} END {if (a>0) b=a*1.2; else b=a*0.8; print b}' )
		fi
		printf "${sample}\t${group}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${ylimmin}\t${ylimmax}\n" >> data/filenames_${ID}.txt
	done < ${samplefile}
	
	if [[ ${htstart} != "" ]] && [[ ${htwidth} != "" ]]; then
		printf "${htstart}\n" | awk -F"," '{for(i=1;i<=NF;i++) print $i}' > data/htstart_${ID}.txt
		printf "${htwidth}\n" | awk -F"," '{for(i=1;i<=NF;i++) print $i}' > data/htwidth_${ID}.txt
		printf "\nPlotting browser on ${ID} with higlights\n\n"
		Rscript --vanilla /grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/R_browser_input_files.r data/filenames_${ID}.txt data/genes_in_locus_${ID}.gff data/tes_in_locus_${ID}.bed ${ID} data/htstart_${ID}.txt data/htwidth_${ID}.txt
	else
		printf "\nPlotting browser on ${ID} without higlights\n\n"
		Rscript --vanilla /grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/R_browser_input_files.r data/filenames_${ID}.txt data/genes_in_locus_${ID}.gff data/tes_in_locus_${ID}.bed ${ID}
	fi
	rm -f data/*${ID}*
done < ${loci}

###################################################

printf "\nScript finished\n"
