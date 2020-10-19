#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 12
#$ -l m_mem_free=12G
#$ -l tmp_free=100G
#$ -o ChIPanalysis.log
#$ -j y
#$ -N ChIPanalysis

usage="
##### Script for Maize code ChIP data analysis
#####
##### Argument #1: samplefile containing the samples to compare and in 5 tab-delimited columns:
##### Line, Tissue, Mark, Rep, PE or SE
##### Argument #2: bedfiles containing the regions that want to be ploted over
##### Safest to use a full path to the region file
##### 
##### It calls broad (for H3K4me2) or narrow (for H3K4me3 and H3K27ac) peaks with Macs2,
##### creates bigwig files (log2 FC vs Input) and makes heatmaps and metaplots with deeptools on the provided regions
##### If the replicates are to be merged, the Rep value in the samplefile should be 'merged'
#####
##### Requirements: samtools, bedtools, deeptools, macs2, idr
"

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS
export samplefile=$1
export regionfile=$2

if [ $# -eq 0 ]; then
	printf "$usage\n"
	exit 1
fi

if [[ "$1" == "help" ]]; then
	printf "$usage\n"
	exit 1
fi

#### To move into ChIP directory if script was launched from MaizeCode folder
path=$(pwd)
curdir=${path##*/}
if [[ $curdir != "ChIP" ]]; then
	cd ChIP
fi

if [ ! -d ./peaks ]; then
	mkdir ./peaks
fi

if [ ! -d ./plots ]; then
	mkdir ./plots
fi


#############################################################################################
########################################### PART1 ###########################################
###################################### Preparing files ######################################
#############################################################################################

linelist=()
tissuelist=()
marklist=()
samplelist=()
bwlist=()
while read line tissue mark rep paired
do
	#### To merge bam files of replicates if chosen and not already exisiting
	if [[ $rep == "merged" ]]; then
		name=${line}_${tissue}_${mark}_${rep}
		input=${line}_${tissue}_Input_${rep}
		tmpname=${line}_${tissue}_${mark}
		tmpinput=${line}_${tissue}_Input
		if [ ! -f mapped/rmdup_${name}.bam ]; then
			printf "\nMerging replicates of $tmpname\n"
			samtools merge -f -@ $threads mapped/temp_${tmpname}.bam mapped/rmdup_${tmpname}_Rep1.bam mapped/rmdup_${tmpname}_Rep2.bam
			samtools sort -@ $threads -o mapped/rmdup_${tmpname}_merged.bam mapped/temp_${tmpname}.bam
			rm -f mapped/temp_${tmpname}.bam
			samtools index -@ $threads mapped/rmdup_${tmpname}_merged.bam
		fi
		if [ ! -f mapped/rmdup_${input}.bam ]; then
			printf "\nMerging replicates of $tmpinput\n"
			samtools merge -@ $threads mapped/temp_${tmpinput}.bam mapped/rmdup_${tmpinput}_Rep1.bam mapped/rmdup_${tmpinput}_Rep2.bam
			samtools sort -@ $threads -o mapped/rmdup_${tmpinput}_merged.bam mapped/temp_${tmpinput}.bam
			rm -f mapped/temp_${tmpinput}.bam
			samtools index -@ $threads mapped/rmdup_${tmpinput}_merged.bam
		fi
	else
		name=${line}_${tissue}_${mark}_${rep}
		input=${line}_${tissue}_Input_${rep}
	fi
	#### To append the lista (used for labels and other stuff)
	samplelist+=("$name")
	linelist+=("$line")
	tissuelist+=("$tissue")
	marklist+=("$mark")
	#### To call either broad or narrow peaks if not disabled and not already exisiting
	case "$mark" in
		H3K4me1) peaktype=broad;;
		H3K4me3) peaktype=narrow;;
		H3K27ac) peaktype=narrow;;
	esac
	if [[ $paired == "PE" ]]; then
		if [[ $peaktype == "broad" ]] && [ ! -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling broad peaks for PE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAMPE -g 2.2e9 -n ${name} --keep-dup "all" --outdir peaks/ --tempdir $TMPDIR --broad
		elif [[ $peaktype == "narrow" ]] && [ ! -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling narrow peaks for PE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAMPE -g 2.2e9 -n ${name} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR
		elif [ -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\n$peaktype peaks already called for $name\n"
		else
			printf "\nSomething is wrong! Check usage by running the script without arguments\n"
			exit 1
		fi
	elif [[ $paired == "SE" ]]; then
		if [[ $peaktype == "broad" ]] && [ ! -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling broad peaks for SE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAM -g 2.2e9 -n ${name} --keep-dup "all" --outdir peaks/ --tempdir $TMPDIR --nomodel --broad
		elif [[ $peaktype == "narrow" ]] && [ ! -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\nCalling narrow peaks for SE $name (vs $input) with macs2 version:\n"
			macs2 --version
			macs2 callpeak -t mapped/rmdup_${name}.bam -c mapped/rmdup_${input}.bam -f BAM -g 2.2e9 -n ${name} --keep-dup "all" --call-summits --outdir peaks/ --tempdir $TMPDIR --nomodel
		elif [[ $peaktype != "narrow" ]] && [[ $peaktype != "broad" ]]; then
			printf "\nType of peaks to call is missing: broad or narrow?\n"
			exit 1
		elif [ -f peaks/${name}_peaks.${peaktype}Peak ]; then
			printf "\n$peaktype peaks already called for $name\n"
		else
			printf "\nSomething is wrong! Check usage by running the script without arguments\n"
			exit 1
		fi
	else
		printf "\nData format missing: paired-end (PE) or single-end (SE)?\n"
		exit 1
	fi
	#### To create bw files if not already exisiting
	if [ ! -f deeptools/${name}.bw ]; then
		printf "\nMaking bigwig files for $name with deeptools version:\n"
		deeptools --version
		bamCompare -b1 mapped/rmdup_${name}.bam -b2 mapped/rmdup_${input}.bam -o deeptools/${name}.bw -p $threads --binSize 1 --scaleFactorsMethod "None" --normalizeUsing CPM
	else
		printf "\nBigwig file for $name already exists! If you want to make a different one (e.g. if you have changed some parameters), delete the existing one or rename it.\n"
	fi
	#### To append the list of bw files to be used by deeptools
	bwlist+=("deeptools/${name}.bw")
done < $samplefile


#############################################################################################
########################################### PART2 ###########################################
################################### Making region files #####################################
#############################################################################################

#### To get IDR analysis on replicates (if they were both present in the samplefile

uniqline=($(printf "%s\n" "${linelist[@]}" | sort -u))
uniqtissue=($(printf "%s\n" "${tissuelist[@]}" | sort -u))
uniqmark=($(printf "%s\n" "${marklist[@]}" | sort -u))

for line in ${uniqline[@]}
do
	for tissue in ${uniqtissue[@]}
	do
		for mark in ${uniqmark[@]}
		do
			case "$mark" in
				H3K4me1) peaktype=broad;;
				H3K4me3) peaktype=narrow;;
				H3K27ac) peaktype=narrow;;
			esac
			if [[ " ${samplelist[@]} " =~ " ${line}_${tissue}_${mark}_Rep1 " ]] && [[ " ${samplelist[@]} " =~ " ${line}_${tissue}_${mark}_Rep2 " ]]
				printf "\nDoing IDR analysis on both replicates from ${line}_${tissue}_${mark} ($peaktype peaks) with idr version:\n"
				idr --version
				idr --input-file-type ${peaktype}Peak --output-file-type ${peaktype}Peak --samples ${line}_${tissue}_${mark}_Rep1.${peaktype}Peak ${line}_${tissue}_${mark}_Rep2.${peaktype}Peak -o peaks/idr_${line}_${tissue}_${mark}.${peaktype}Peak -l reports/idr_${line}_${tissue}_${mark}.log --plot --use-best-multisummit-IDR
			fi
		done
	done
done

# #### To compare peaks among the different marks and make groups based on how many marks are shared

# #### The name of the files is a combination of the samplefile name and the regionfile name

# tmp1=${samplefile##*/}
# tmp2=${tmp1%.*}
# tmp3=${regionfile##*/}
# tmp4=${tmp3%.*}
# analysisname="${tmp2}_on_${tmp4}"

# for name in ${samplelist[@]}
# do
	# awk -v OFS="\t" -v n=$name '{print $1,$2-1,$3,s"_"NR}' peaks/${name}_peaks.xls > peaks/peaks_${name}.bed
# done
# cat peaks/peaks_*.bed | sort -k1,1n -k2,2n > peaks/peaks_${analysisname}.bed

# bedtools merge -i peaks/peaks_${analysisname}.bed -c 4 -o distinct > peaks/peaks_${analysisname}.bed

# #### To get the distance between the peaks in the different groups


# #############################################################################################
# ########################################### PART3 ###########################################
# ####################################### Making plots ########################################
# #############################################################################################


# #### To make heatmaps and profiles with deeptools
# #### By default, it does both scale-regions and reference-point on start of bedfile provided
# #### By default, it does heatmap on all the data, heatmap with 5 kmeans, and corresponding profiles
# #### Probably need to edit many parameters depending on the purpose of the analysis

# #### The name of the files and plots is a combination of the samplefile name and the regionfile name

# printf "\nDoing analysis for $analysisname with deeptools version:\n"
# deeptools --version

# #### Computing the matrix
# if [ ! -f deeptools/regions_${analysisname}.gz ]; then
	# printf "\nComputing scale-regions matrix for $analysisname\n"
	# computeMatrix scale-regions -R $regionfile -S ${bwlist[@]} -bs 50 -b 2000 -a 2000 -m 5000 -p $threads -o deeptools/regions_${analysisname}.gz
# elif [ ! -f deeptools/tss_${analysisname}.gz ]; then
	# printf "\nComputing reference-point on TSS matrix for $analysisname\n"
	# computeMatrix reference-point --referencePoint "TSS" -R $regionfile -S ${bwlist[@]} -bs 50 -b 2000 -a 6000 -p $threads -o deeptools/tss_${analysisname}.gz
# else
	# printf "\nMatrix already computed for ${analysisname}! If you want to make a different one (e.g. if you have changed some parameters), delete the existing one or rename it.\n"
# fi

# #### Ploting heatmaps
# printf "\nPlotting full heatmap for scale-regions of $analysisname\n"
# plotHeatmap -m deeptools/regions_${analysisname}.gz -out plots/${analysisname}_heatmap_regions.pdf --sortRegions descend --sortUsing mean --samplesLabel ${samplelist[@]} --colorMap 'seismic'
# printf "\nPlotting heatmap for scale-regions of $analysisname split in 5 kmeans\n"
# plotHeatmap -m deeptools/regions_${analysisname}.gz -out plots/${analysisname}_heatmap_regions_k5.pdf --sortRegions descend --sortUsing mean --samplesLabel ${samplelist[@]} --colorMap 'seismic' --kmeans 5 --outFileSortedRegions deeptools/${analysisname}_sortedregions_k5.txt

# printf "\nPlotting full heatmap for reference-point TSS of $analysisname\n"
# plotHeatmap -m deeptools/tss_${analysisname}.gz -out plots/${analysisname}_heatmap_tss.pdf --sortRegions descend --sortUsing region_length --samplesLabel ${samplelist[@]} --colorMap 'seismic'
# printf "\nPlotting heatmap for reference-point TSS of $analysisname split in 5 kmeans\n"
# plotHeatmap -m deeptools/tss_${analysisname}.gz -out plots/${analysisname}_heatmap_tss_k5.pdf --sortRegions descend --sortUsing region_length --samplesLabel ${samplelist[@]} --colorMap 'seismic' --kmeans 5 --outFileSortedRegions deeptools/${analysisname}_sortedtss_k5.txt

# #### Plotting Metaplot profiles
# printf "\nPlotting metaplot profiles for scale-regions of $analysisname\n"
# plotProfile -m deeptools/regions_${analysisname}.gz -out plots/${analysisname}_profiles_regions.pdf --plotType lines --averageType mean --perGroup
# printf "\nPlotting metaplot profiles for scale-regions of $analysisname split in 5 kmeans\n"
# plotProfile -m deeptools/regions_${analysisname}.gz -out plots/${analysisname}_profiles_regions_k5.pdf --plotType lines --averageType mean --perGroup --kmeans 5

# printf "\nPlotting metaplot profiles for reference-point TSS of $analysisname\n"
# plotProfile -m deeptools/tss_${analysisname}.gz -out plots/${analysisname}_profiles_tss.pdf --plotType lines --averageType mean --perGroup
# printf "\nPlotting metaplot profiles for reference-point TSS of $analysisname split in 5 kmeans\n"
# plotProfile -m deeptools/tss_${analysisname}.gz -out plots/${analysisname}_profiles_tss_k5.pdf --plotType lines --averageType mean --perGroup --kmeans 5

# #### When done this way, the 5 kmeans regions in heatmap and profiles are not going to be the same. 
# #### To have the same regions, make a new matrix using the region file coming from the --outFileSortedRegions (e.g. deeptools/${analysisname}_sortedtss_k5.txt)
# #### You can then keep its order (--sortUsing keep) if required


#############################################################################################
########################################### MISC ############################################
#############################################################################################

###### To make a test samplefile for (B73 endosperm H3K4me1)

# printf "B73\tendosperm\tH3K4me1\tRep1\tPE\nB73\tendosperm\tH3K4me1\tRep2\tPE\n" > test_analysis_samplefile.txt

###### To make the samplefile for B73_endosperm

# printf "B73\tendosperm\tH3K4me1\tRep1\tPE\nB73\tendosperm\tH3K4me1\tRep2\tPE\nB73\tendosperm\tH3K4me1\tmerged\tPE\nB73\tendosperm\tH3K4me3\tRep1\tPE\nB73\tendosperm\tH3K4me3\tRep2\tPE\nB73\tendosperm\tH3K4me3\tmerged\tPE\nB73\tendosperm\tH3K4ac\tRep1\tPE\nB73\tendosperm\tH3K4ac\tRep2\tPE\nB73\tendosperm\tH3K4ac\tmerged\tPE\n" > B73_endosperm_analysis_samplefile.txt

###### To create a regionfile containing several groups of regions (e.g. gene_list1, gene_list2 and gene_list3)

# printf "#Name_gene_list1" > regionfile.bed
# cat gene_list1.bed >> regionfile.bed
# printf "#Name_gene_list2" >> regionfile.bed
# cat gene_list2.bed >> regionfile.bed
# printf "#Name_gene_list3" >> regionfile.bed
# cat gene_list3.bed >> regionfile.bed

#############################################################################################


printf "\nScript finished successfully!\n"		