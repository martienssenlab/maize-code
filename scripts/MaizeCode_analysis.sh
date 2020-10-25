#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 8
#$ -l m_mem_free=10G
#$ -l tmp_free=50G
#$ -o maizecodeanalysis.log
#$ -j y
#$ -N maizecodeanalysis

usage="
##### Script for Maize code data analysis
#####
##### sh MaiCode_analysis.sh -f samplefile -r regionfile [-s]
#####	-f: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### 		Line, Tissue, Sample, Rep (Rep1, Rep2 or merged), PE or SE
##### 	-r: bedfiles containing the regions that want to be ploted over
##### 		(safest to use a full path to the region file)
#####		If no region file is given, the analysis will behave as if -s was set
#####	-s: If set, the script stops after calling peaks for ChIP samples or DEG for RNA samples and making bigwigs
#####		Set it last to avoid problems
##### 	-h: help, returns usage
##### 
##### It sends each type of sample to its specific analysis file (MaizeCode_ChIP_analysis.sh or MaizeCode_RNA_analysis.sh)
##### Then combines into deeper analysis (to detail...)
#####
##### Requirements: samtools, bedtools, deeptools, macs2, idr, R (+R packages: ggplot2,readr,UpSetR)
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

while getopts ":f:r:sh" opt; do
	case $opt in
		h) 	printf "$usage\n"
			exit 0;;
		f) 	export samplefile=${OPTARG};;
		r)	export regionfile=${OPTARG};;
		s)	printf "\nOption to stop after bigwigs selected\n"
			export keepgoing="STOP";;
		*)	printf "$usage\n"
			exit 1;;
	esac
done
shift $((OPTIND - 1))

if [ ! $samplefile ]; then
	printf "Missing samplefile\n"
	printf "$usage\n"
	exit 1
fi

if [ ! $regionfile ]; then
	printf "Regionfile missing or not bed format.\nAnalysis will be stopped without deeper analysis\nIf that was not intended, check usage.\n"
	keepgoing="STOP"
fi

#############################################################################################
########################################### PART1 ###########################################
######################## Send each sample to specific analysis type #########################
#############################################################################################

tmp1=${samplefile##*/}
samplename=${tmp1%_analysis*}

if [ -f temp_${samplename}_ChIP.txt ]; then
	rm temp_${samplename}_ChIP.txt
fi

if [ -f temp_${samplename}_RNA.txt ]; then
	rm temp_${samplename}_RNA.txt
fi

datatype_list=()
while read line tissue sample rep paired
do
	case "$sample" in
		H*|Input) datatype="ChIP";;
		*RNA*|RAMPAGE) datatype="RNA";;
	esac
	if [[ "$datatype" == "ChIP" ]]; then
		printf "$line\t$tissue\t$sample\t$rep\t$paired\n" >> $datatype/temp_${samplename}_ChIP.txt
	elif [[ "$datatype" == "RNA" ]]; then
		printf "$line\t$tissue\t$sample\t$rep\t$paired\n" >> $datatype/temp_${samplename}_RNA.txt
	else
		printf "Type of data unknown!\n"
		exit 1
	fi
	datatype_list+=("$datatype")
done < $samplefile

uniq_datatype_list=($(printf "%s\n" "${datatype_list[@]}" | sort -u))

pids=()
for datatype in ${uniq_datatype_list[@]}
do
	printf "\nRunning $datatype analysis script\n"
	cd $datatype
	qsub -N ${datatype} -o ${datatype}.log ~/data/Scripts/MaizeCode_${datatype}_analysis.sh -f temp_${samplename}_${datatype}.txt -r ${regionfile} &
	pids+=("$!")
	cd ..
done

#### Wait for the other scripts to finish
printf "\nWaiting for the datatype analysis scripts to finish\n"
wait $pids || { printf "\nThere was an error during analysis\n" >&2; exit 1; }

if [[ $keepgoing == "STOP" ]]; then
	printf "\nScript finished successfully without deeper analysis\n"		
	exit 0
fi	

#############################################################################################
########################################### PART2 ###########################################
############################### Combining ChIP and RNA data #################################
#############################################################################################





#############################################################################################
########################################### PART3 ###########################################
####################################### Making plots ########################################
#############################################################################################


#### To make heatmaps and profiles with deeptools
#### By default, it does both scale-regions and reference-point on start of bedfile provided
#### By default, it does heatmap on all the data, heatmap with 5 kmeans, and corresponding profiles
#### Probably need to edit many parameters depending on the purpose of the analysis

# printf "\nDoing analysis for $analysisname with deeptools version:\n"
# deeptools --version

# #### Computing the matrix
# if [ ! -f tracks/regions_${analysisname}.gz ]; then
	# printf "\nComputing scale-regions matrix for $analysisname\n"
	# computeMatrix scale-regions -R $regionfile -S ${bw_list[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o tracks/regions_${analysisname}.gz
# fi
# if [ ! -f tracks/tss_${analysisname}.gz ]; then
	# printf "\nComputing reference-point on TSS matrix for $analysisname\n"
	# computeMatrix reference-point --referencePoint "TSS" -R $regionfile -S ${bw_list[@]} -bs 50 -b 2000 -a 6000 -p $threads -o tracks/tss_${analysisname}.gz
# fi

# #### Ploting heatmaps
# printf "\nPlotting full heatmap for scale-regions of $analysisname\n"
# plotHeatmap -m tracks/regions_${analysisname}.gz -out plots/${analysisname}_heatmap_regions.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sample_list[@]} --colorMap 'seismic'
# printf "\nPlotting heatmap for scale-regions of $analysisname split in 3 kmeans\n"
# plotHeatmap -m tracks/regions_${analysisname}.gz -out plots/${analysisname}_heatmap_regions_k3.pdf --sortRegions descend --sortUsing mean --samplesLabel ${sample_list[@]} --colorMap 'seismic' --kmeans 3 --outFileSortedRegions tracks/${analysisname}_sortedregions_k3.txt

# printf "\nPlotting full heatmap for reference-point TSS of $analysisname\n"
# plotHeatmap -m tracks/tss_${analysisname}.gz -out plots/${analysisname}_heatmap_tss.pdf --sortRegions descend --sortUsing region_length --samplesLabel ${sample_list[@]} --colorMap 'seismic'
# printf "\nPlotting heatmap for reference-point TSS of $analysisname split in 3 kmeans\n"
# plotHeatmap -m tracks/tss_${analysisname}.gz -out plots/${analysisname}_heatmap_tss_k3.pdf --sortRegions descend --sortUsing region_length --samplesLabel ${sample_list[@]} --colorMap 'seismic' --kmeans 3 --outFileSortedRegions tracks/${analysisname}_sortedtss_k3.txt

# #### Plotting Metaplot profiles
# printf "\nPlotting metaplot profiles for scale-regions of $analysisname\n"
# plotProfile -m tracks/regions_${analysisname}.gz -out plots/${analysisname}_profiles_regions.pdf --plotType lines --averageType mean --perGroup --samplesLabel ${sample_list[@]}
# printf "\nPlotting metaplot profiles for scale-regions of $analysisname split in 5 kmeans\n"
# plotProfile -m tracks/regions_${analysisname}.gz -out plots/${analysisname}_profiles_regions_k5.pdf --plotType lines --averageType mean --perGroup --samplesLabel ${sample_list[@]} --kmeans 5

# printf "\nPlotting metaplot profiles for reference-point TSS of $analysisname\n"
# plotProfile -m tracks/tss_${analysisname}.gz -out plots/${analysisname}_profiles_tss.pdf --plotType lines --averageType mean --perGroup --samplesLabel ${sample_list[@]}
# printf "\nPlotting metaplot profiles for reference-point TSS of $analysisname split in 5 kmeans\n"
# plotProfile -m tracks/tss_${analysisname}.gz -out plots/${analysisname}_profiles_tss_k5.pdf --plotType lines --averageType mean --perGroup --samplesLabel ${sample_list[@]} --kmeans 5

#### When done this way, the 5 kmeans regions in heatmap and profiles are not going to be the same. 
#### To have the same regions, make a new matrix using the region file coming from the --outFileSortedRegions (e.g. tracks/${analysisname}_sortedtss_k5.txt)
#### You can then keep its order (--sortUsing keep) if required


#############################################################################################
########################################### MISC ############################################
#############################################################################################

###### To make a test samplefile for (B73 endosperm H3K4me1)

# printf "B73\tendosperm\tH3K4me1\tRep1\tPE\nB73\tendosperm\tH3K4me1\tRep2\tPE\n" > test_analysis_samplefile.txt

###### To make the samplefile for B73_endosperm

# printf "B73\tendosperm\tH3K4me1\tRep1\tPE\nB73\tendosperm\tH3K4me1\tRep2\tPE\nB73\tendosperm\tH3K4me1\tmerged\tPE\nB73\tendosperm\tH3K4me3\tRep1\tPE\nB73\tendosperm\tH3K4me3\tRep2\tPE\nB73\tendosperm\tH3K4me3\tmerged\tPE\nB73\tendosperm\tH3K27ac\tRep1\tPE\nB73\tendosperm\tH3K27ac\tRep2\tPE\nB73\tendosperm\tH3K27ac\tmerged\tPE\n" > B73_endosperm_analysis_samplefile.txt

###### To create a regionfile containing several groups of regions (e.g. gene_list1, gene_list2 and gene_list3)

# printf "#Name_gene_list1" > regionfile.bed
# cat gene_list1.bed >> regionfile.bed
# printf "#Name_gene_list2" >> regionfile.bed
# cat gene_list2.bed >> regionfile.bed
# printf "#Name_gene_list3" >> regionfile.bed
# cat gene_list3.bed >> regionfile.bed

#############################################################################################


printf "\nScript finished successfully!\n"
