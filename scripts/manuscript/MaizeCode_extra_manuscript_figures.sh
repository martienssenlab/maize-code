#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 24
#$ -l m_mem_free=4G
#$ -l tmp_free=8G
#$ -o final.log
#$ -j y
#$ -N final

set -e -o pipefail

export threads=$NSLOTS

printf "\n\n"
date
printf "\n"

################################################################
### To plot the heatmap and profile of B73 on all genes (Fig 1a)
################################################################

# line="B73"
# ref="B73_v5"
# matrix="/grid/martienssen/home/jcahn/norepl/projects/MNase/combined/matrix/all_genes_regions_tot_B73_on_B73_v5_all_genes_for_H3K27ac.gz"
# nameout="Fig1a_B73_ChIP_on_all_genes"
# matrixout="manuscript/matrix/${nameout}.gz"

# samples=()
# labels=()
# colorsmap=()
# colorsline=()
# for tissue in ears cn roots endosperm
# do
	# samples+=("${line}_${tissue}_H3K27ac_merged" "${line}_${tissue}_H3K4me1_merged" "${line}_${tissue}_H3K4me3_merged" "${line}_${tissue}_RNAseq_merged_plus" "${line}_${tissue}_RAMPAGE_merged_plus" "${line}_${tissue}_MNase_merged")
	# labels+=("${tissue}_H3K27ac" "${tissue}_H3K4me1" "${tissue}_H3K4me3" "${tissue}_RNAseq" "${tissue}_RAMPAGE" "${tissue}_MNase")
	# colorsmap+=("Reds" "Blues" "YlOrBr" "Purples" "Greys" "Greens")
	# colorsline+=("firebrick" "cornflowerblue" "goldenrod" "darkmagenta" "dimgray" "darkgreen")
# done

# printf "subsetting matrix for ${nameout}\n"
# computeMatrixOperations subset -m ${matrix} -o ${matrixout} --samples ${samples[@]}
# printf "\nGetting Z values for ${nameout}\n"
# computeMatrixOperations dataRange -m ${matrixout} > manuscript/matrix/values_z_${nameout}.txt
# printf "\nGetting scales for ${nameout}\n"
# zmins=()
# zmaxs=()
# for sample in ${samples[@]}
# do
	# mini=$(grep ${sample} manuscript/matrix/values_z_${nameout}.txt | awk '{print $5}')
	# maxi=$(grep ${sample} manuscript/matrix/values_z_${nameout}.txt | awk '{print $6}')
	# test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
	# if [[ ${test} == "yes" ]]; then
		# zmins+=("0")
		# zmaxs+=("0.005")
	# else
		# zmins+=("${mini}")
		# zmaxs+=("${maxi}")
	# fi
# done
# printf "\nGetting Y values for ${nameout}\n"
# plotProfile -m ${matrixout} -out manuscript/plots/temp_${nameout}_profile.pdf --averageType mean --outFileNameData manuscript/matrix/values_y_${nameout}.txt
# rm -f manuscript/plots/temp_${nameout}_profile.pdf
# printf "\nGetting Y scales for ${nameout}\n"
# ymins=()
# ymaxs=()
# for sample in ${samples[@]}
# do
	# ymini=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
	# ymaxi=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
	# test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
	# if [[ ${test} == "yes" ]]; then
		# ymins+=("0")
		# ymaxs+=("0.01")
	# else
		# ymins+=("${ymini}")
		# ymaxs+=("${ymaxi}")
	# fi
# done
# printf "\nPlotting heatmap for ${nameout}\n"
# plotHeatmap -m ${matrixout} -out manuscript/plots/${nameout}_heatmap.pdf --sortRegions descend --sortUsing mean --samplesLabel ${labels[@]} --regionsLabel "" --colorMap ${colorsmap[@]} --zMin ${zmins[@]} --zMax ${zmaxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --whatToShow "heatmap only" --boxAroundHeatmaps "no"
# printf "\nPlotting heatmap for ${nameout}\n"
# plotHeatmap -m ${matrixout} -out manuscript/plots/${nameout}_heatmap_full.pdf --sortRegions descend --sortUsing mean --samplesLabel ${labels[@]} --regionsLabel "" --colorMap ${colorsmap[@]} --zMin ${zmins[@]} --zMax ${zmaxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --boxAroundHeatmaps "no"
# printf "\nPlotting profile for ${nameout}\n"
# plotProfile -m ${matrixout} -out manuscript/plots/${nameout}_profile.pdf --samplesLabel ${labels[@]} --averageType mean --colors ${colorsline[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --numPlotsPerRow ${#samples[@]}


########################################
### To plot the heatmap on OCRs (Fig 1b)
########################################

# pathv4="/grid/martienssen/home/jcahn/norepl/projects/maizecode_v4"
# colorsmap=("Reds" "Blues" "YlOrBr")
# colorsline+=("firebrick" "cornflowerblue" "goldenrod")
# bw_list=("${pathv4}/ChIP/tracks/B73_ears_H3K27ac_merged.bw" "${pathv4}/ChIP/tracks/B73_ears_H3K4me1_merged.bw" "${pathv4}/ChIP/tracks/B73_ears_H3K4me3_merged.bw")
# label_list=("H3K27ac" "H3K4me1" "H3K4me3")

# awk -v OFS="\t" 'NR>1 && $1~/^[1-9]/ {print $1,$2,$3,$8}' manuscript/Sun_2020_S2_ear_OCRs.txt | sort -k1,1n -k2,2n > manuscript/matrix/OCRs_ears.bed

# nameout="Fig1b_ears_OCRs"
# region_label=()
# for type in LoOCR dOCR
# do
	# awk -v t=${type} '$4~t' manuscript/matrix/OCRs_ears.bed > manuscript/matrix/OCRs_ears_${type}.bed
	# n=$(wc -l manuscript/matrix/OCRs_ears_${type}.bed | awk '{print $1}')
	# region_label+=("${type}(${n})")
# done
# printf "\nComputing matrix for ${nameout}\n"
# computeMatrix scale-regions -q --missingDataAsZero --skipZeros -R manuscript/matrix/OCRs_ears_LoOCR.bed manuscript/matrix/OCRs_ears_dOCR.bed -S ${bw_list[@]} -bs 10 -b 2000 -a 2000 -m 1000 -p ${threads} -o manuscript/matrix/${nameout}.gz
# printf "\nGetting y values for ${nameout}\n"
# plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_temp_profile.pdf --samplesLabel ${label_list[@]} --averageType mean --outFileNameData manuscript/matrix/values_y_${nameout}.txt
# rm -f manuscript/plots/${nameout}_temp_profile.pdf
# printf "\nGetting y scales for ${nameout}\n"
# ymins=()
# ymaxs=()
# for sample in ${label_list[@]}
# do
	# ymini=$(grep $sample manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*2; else a=m*0.5; print a}')
	# ymaxi=$(grep $sample manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m*2}')
	# ymins+=("$ymini")
	# ymaxs+=("$ymaxi")
# done
# printf "\nGetting z values for ${nameout}\n"
# computeMatrixOperations dataRange -m manuscript/matrix/${nameout}.gz > manuscript/matrix/values_z_${nameout}.txt
# printf "\nGetting z scales for ${nameout}\n"
# zmins=()
# zmaxs=()
# totnb=${#label_list[*]}
# for (( i=1; i<=${totnb}; i++ ))
# do 
	# min1=$(awk -v i=$i 'NR==(i+1) {print $5}' manuscript/matrix/values_z_${nameout}.txt)
	# max1=$(awk -v i=$i 'NR==(i+1) {print $6}' manuscript/matrix/values_z_${nameout}.txt)
	# zmins+=("$min1")
	# zmaxs+=("$max1")
# done
# printf "\nPlotting heatmap for ${nameout}\n"
# plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_heatmap.pdf --sortRegions descend --sortUsing mean --sortUsingSamples 1 --samplesLabel ${label_list[@]} --regionsLabel ${region_label[@]} --colorMap ${colorsmap[@]} --interpolationMethod 'bilinear' --yMin ${ymins[@]} --yMax ${ymaxs[@]} --zMin ${zmins[@]} --zMax ${zmaxs[@]} --whatToShow "heatmap only" --boxAroundHeatmaps "no"
# plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_heatmap_full.pdf --sortRegions descend --sortUsing mean --sortUsingSamples 1 --samplesLabel ${label_list[@]} --regionsLabel ${region_label[@]} --colorMap ${colorsmap[@]} --interpolationMethod 'bilinear' --yMin ${ymins[@]} --yMax ${ymaxs[@]} --zMin ${zmins[@]} --zMax ${zmaxs[@]} --boxAroundHeatmaps "no"
# printf "\nPlotting profile for ${nameout}\n"
# plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_profile.pdf --samplesLabel ${label_list[@]} --averageType mean --yMin ${ymins[@]} --yMax ${ymaxs[@]} --numPlotsPerRow ${#label_list[@]}

###########################################################################
### To plot the heatmap of intersection between enhancers and OCRs (Fig 1c)
###########################################################################

# pathv4="/grid/martienssen/home/jcahn/norepl/projects/maizecode_v4"
# colorsmap=("Reds" "Blues" "YlOrBr")
# colorsline+=("firebrick" "cornflowerblue" "goldenrod")
# bw_list=("${pathv4}/ChIP/tracks/B73_ears_H3K27ac_merged.bw" "${pathv4}/ChIP/tracks/B73_ears_H3K4me1_merged.bw" "${pathv4}/ChIP/tracks/B73_ears_H3K4me3_merged.bw")
# label_list=("H3K27ac" "H3K4me1" "H3K4me3")

# awk -v OFS="\t" 'NR>1 && $1~/^[1-9]/ {print $1,$2,$3,$8}' manuscript/Sun_2020_S2_ear_OCRs.txt | sort -k1,1n -k2,2n > manuscript/matrix/OCRs_ears.bed

# nameout="Fig1c_ears_OCRs_intersect"
# rm -f manuscript/matrix/B73_v4_ears_enhancers.bed
# rm -f manuscript/matrix/${nameout}_overlap.txt
# rm -f manuscript/matrix/${nameout}_no_overlap.txt
# for type in genic promoter terminator distal_upstream distal_downstream
# do
	# awk -v OFS="\t" -v t=${type} 'NR>1 {print $1,$2,$3,t}' ${pathv4}/combined/peaks/complete_enhancers_${type}_B73_ears_B73_v4_on_B73_v4_all_genes.txt >> manuscript/matrix/B73_v4_ears_enhancers.bed
# done
# bedtools intersect -wao -a manuscript/matrix/B73_v4_ears_enhancers.bed -b manuscript/matrix/OCRs_ears.bed | awk -v OFS="\t" -v x=${nameout} '{if ($9!=0) printf $1"\t"$2"\t"$3"\n" >> "manuscript/matrix/"x"_overlap.txt"; else printf $1"\t"$2"\t"$3"\n" >> "manuscript/matrix/"x"_no_overlap.txt" }'
# nbyes=$(wc -l manuscript/matrix/${nameout}_overlap.txt | awk '{print $1}')
# nbno=$(wc -l manuscript/matrix/${nameout}_no_overlap.txt | awk '{print $1}')
# regionlabel="Overlap_OCRs(${nbyes}) No_overlap(${nbno})"
# printf "\nComputing matrix for ${nameout}\n"
# computeMatrix scale-regions -q --missingDataAsZero --skipZeros -R manuscript/matrix/${nameout}_overlap.txt manuscript/matrix/${nameout}_no_overlap.txt -S ${bw_list[@]} -bs 10 -b 2000 -a 2000 -m 1000 -p ${threads} -o manuscript/matrix/${nameout}.gz
# printf "\nGetting y values for ${nameout}\n"
# plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_temp_profile.pdf --samplesLabel ${label_list[@]} --averageType mean --outFileNameData manuscript/matrix/values_y_${nameout}.txt
# rm -f manuscript/plots/${nameout}_temp_profile.pdf
# printf "\nGetting y scales for ${nameout}\n"
# ymins=()
# ymaxs=()
# for sample in ${label_list[@]}
# do
	# ymini=$(grep $sample manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*2; else a=m*0.5; print a}')
	# ymaxi=$(grep $sample manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {print m*2}')
	# ymins+=("$ymini")
	# ymaxs+=("$ymaxi")
# done
# printf "\nGetting z values for ${nameout}\n"
# computeMatrixOperations dataRange -m manuscript/matrix/${nameout}.gz > manuscript/matrix/values_z_${nameout}.txt
# printf "\nGetting z scales for ${nameout}\n"
# zmins=()
# zmaxs=()
# totnb=${#label_list[*]}
# for (( i=1; i<=${totnb}; i++ ))
# do 
	# min1=$(awk -v i=$i 'NR==(i+1) {print $5}' manuscript/matrix/values_z_${nameout}.txt)
	# max1=$(awk -v i=$i 'NR==(i+1) {print $6}' manuscript/matrix/values_z_${nameout}.txt)
	# zmins+=("$min1")
	# zmaxs+=("$max1")
# done
# printf "\nPlotting heatmap for ${nameout}\n"
# plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_heatmap.pdf --sortRegions descend --sortUsing mean --sortUsingSamples 1 --samplesLabel ${label_list[@]} --regionsLabel ${regionlabel} --colorMap ${colorsmap[@]} --interpolationMethod 'bilinear' --yMin ${ymins[@]} --yMax ${ymaxs[@]} --zMin ${zmins[@]} --zMax ${zmaxs[@]} --whatToShow "heatmap only" --boxAroundHeatmaps "no"
# plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_heatmap_full.pdf --sortRegions descend --sortUsing mean --sortUsingSamples 1 --samplesLabel ${label_list[@]} --regionsLabel ${regionlabel} --colorMap ${colorsmap[@]} --interpolationMethod 'bilinear' --yMin ${ymins[@]} --yMax ${ymaxs[@]} --zMin ${zmins[@]} --zMax ${zmaxs[@]} --boxAroundHeatmaps "no"
# printf "\nPlotting profile for ${nameout}\n"
# plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_profile.pdf --samplesLabel ${label_list[@]} --averageType mean --yMin ${ymins[@]} --yMax ${ymaxs[@]} --numPlotsPerRow ${#label_list[@]} --regionsLabel ${regionlabel}

#################################### 
### To do the TF Upset plot (Fig 2a)
#################################### 

# line="B73"
# ref="B73_v5"
# analysisname="B73_on_B73_v5_all_genes"
# nameout="SupFig1a_TF_Upset_mix"
# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/B73_v5"
# regionfile="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/shRNA/tracks/B73_v5_all_genes.bed"
# tf_sample_list="GT1 HDZIV6 TU1A TB1 FEA4 KN1"

# ### Using idr file for TB1 and selected peaks for the others

# printf "\nPreparing merged TF peaks file for ${analysisname}\n"
# if [ -s manuscript/data/tmp_peaks_${nameout}.bed ]; then
	# rm -f manuscript/data/tmp_peaks_${nameout}.bed
# fi
# for sample in ${tf_sample_list[@]} H3K27ac
# do
	# case "${sample}" in
		# H3K27ac)	file="combined/peaks/merged_peaks_H3K27ac_${analysisname}.bed";;
		# TB1)	file="TF/peaks/idr_${line}_${sample}.narrowPeak";;
		# *)	file="TF/peaks/selected_peaks_${line}_${sample}.narrowPeak";;
	# esac
	# awk -v OFS="\t" -v s=${sample} '($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/ ) {print $1,$2,$3,s}' ${file} | sort -k1,1 -k2,2n -u >> manuscript/data/tmp_${nameout}.bed
# done

# sort -k1,1 -k2,2n manuscript/data/tmp_${nameout}.bed > manuscript/data/tmp2_${nameout}.bed
# bedtools merge -i manuscript/data/tmp2_${nameout}.bed -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" ' $4 != "H3K27ac" {print $1,$2,$3,"Peak_"NR,$4}'> manuscript/data/tmp3_${nameout}.bed
# ## To get distance to closest gene (and the gene model name)
# printf "\nGetting closest region of ${nameout}\n"
# bedtools closest -a manuscript/data/tmp3_${nameout}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > manuscript/data/${nameout}.bed
# rm -f manuscript/data/tmp*_${nameout}.bed
# ## To create a matrix of peak presence in each sample
# printf "\nCreating matrix file for ${nameout}\n"
# for sample in ${tf_sample_list[@]} H3K27ac
# do
	# printf "${sample}\n" > manuscript/data/temp_col_${nameout}_${sample}.txt
	# awk -v OFS="\t" -v s=${sample} '{if ($0 ~ s) print "1"; else print "0"}' manuscript/data/${nameout}.bed >> manuscript/data/temp_col_${nameout}_${sample}.txt
# done

# ## To group peaks based on their distance (gene body (x=0), promoter (0<x<2kb upstream), terminator (0<x<2kb downstream), distal)
# awk -v OFS="\t" 'BEGIN {printf "PeakID\tDistance\tGroup\n"} {if ($5<-2000) {d="Distal_downstream"; a=-$5} else if ($5<0) {d="Terminator"; a=-$5} else if ($5==0) {d="Gene_body"; a=$5} else if ($5>2000) {d="Distal_upstream"; a=$5} else {d="Promoter"; a=$5} print $4,a,d}' manuscript/data/${nameout}.bed > manuscript/data/temp_col_${nameout}_AAA.txt
# paste manuscript/data/temp_col_${nameout}_*.txt | uniq > manuscript/data/matrix_${nameout}.txt
# rm -f manuscript/data/temp_col_${nameout}_*.txt
# ## To make an Upset plot highlighting peaks in gene bodies
# printf "\nCreating Upset plot for ${nameout} with R version:\n"
# R --version
# Rscript --vanilla MaizeCode_R_Upset_TF.r ${nameout} manuscript/data/matrix_${nameout}.txt

####################################################################################################################
#### To make browser shots at TF loci (Fig 2c)
###################################################################################################################

# printf "tb1\tteosinte_branched1\tZm00001eb054440\n" > manuscript/data/tf2_genes_B73.txt
# printf "fea4\tfasciated_ear4\tZm00001eb280500\n" >> manuscript/data/tf2_genes_B73.txt
# printf "ra1\tramosa1\tZm00001eb312340\n" >> manuscript/data/tf2_genes_B73.txt
# printf "gt1\tgrassy_tillers1\tZm00001eb007950\n" >> manuscript/data/tf2_genes_B73.txt
# printf "tga1\tteosinte_glume_architecture1\tZm00001eb175150\n" >> manuscript/data/tf2_genes_B73.txt

# # for samp in all_tf subset_tf
# for samp in subset_tf
# do
	# nameout="browser_${samp}"
	# samplefile="${samp}_analysis_samplefile.txt"
	# rm -f browser/data/${nameout}_samplefile.txt
	
	# while read data line tissue sample paired ref_dir
	# do
		# if [[ ${data} = RNAseq ]] || [[ ${data} = RAMPAGE ]] || [[ ${data} = shortRNA ]] || [[ ${data} = 24RNA ]]; then
			# case "${data}" in
				# RNAseq)		folder="RNA"	
							# backcolor="#7716ab"
							# trackcolor="#7716ab"
							# fillcolor1="#9067a6"
							# fillcolor2="#9067a6";;
				# RAMPAGE)	folder="RNA"	
							# backcolor="grey10"
							# trackcolor="grey10"
							# fillcolor1="grey50"
							# fillcolor2="grey50";;
				# shortRNA)	folder="shRNA"
							# backcolor="#C45100"
							# trackcolor="#C45100"
							# fillcolor1="#FB9E00"
							# fillcolor2="#FB9E00";;
				# 24RNA)	folder="shRNA"
						# backcolor="#C45100"
						# trackcolor="#C45100"
						# fillcolor1="#FB9E00"
						# fillcolor2="#FB9E00";;
			# esac
			# pathtobw="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/${folder}/tracks"
			# for strand in plus minus
			# do
				# case "${strand}" in
					# plus)	invert="no";;
					# minus)	invert="invert";;
				# esac
				# bw="${pathtobw}/${line}_${tissue}_${sample}_merged_${strand}.bw"
				# printf "${line}_${tissue}\t${sample}_${strand}\t${bw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\n" >> browser/data/${nameout}_samplefile.txt
			# done
		# elif [[ ${data} = ChIP* ]]; then
			# case "${sample}" in 
				# H3K27ac)	pathtobw="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/ChIP/tracks"	
						# backcolor="#EE616E"
						# trackcolor="#EE616E"
						# fillcolor1="#ebd8da"
						# fillcolor2="#ed939b"
						# lims="-0.2;1.1";;
				# H3K4me1)	pathtobw="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/ChIP/tracks"	
						# backcolor="#8D9BEE"
						# trackcolor="#8D9BEE"
						# fillcolor1="#d5d8e8"
						# fillcolor2="#afb8ed"
						# lims="";;
				# H3K4me3)	pathtobw="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/ChIP/tracks"	
						# backcolor="#F1C062"
						# trackcolor="#F1C062"
						# fillcolor1="#e8e1d3"
						# fillcolor2="#ebcd94"
						# lims="";;
				# MNase)	pathtobw="/grid/martienssen/home/jcahn/norepl/projects/MNase/ChIP/tracks"
						# backcolor="grey"
						# trackcolor="grey"
						# fillcolor1="grey"
						# fillcolor2="grey"
						# lims="";;
			# esac
			# bw="${pathtobw}/${line}_${tissue}_${sample}_merged.bw"
			# printf "${line}_${tissue}\t${sample}\t${bw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\tno\t${lims}\n" >> browser/data/${nameout}_samplefile.txt
		# elif [[ ${data} = TF* ]]; then
			# # pathtobw="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/TF/tracks"
			# pathtobw="/grid/martienssen/home/jcahn/norepl/projects/MNase/TF/tracks"
			# name=${data##TF_}
			# case "${name}" in 
				# TB1)	backcolor="red"
						# trackcolor="red"
						# fillcolor1="red"
						# fillcolor2="red"
						# lims="-2.5;4";;
						# # lims="-2.4;5.5";;
						# # lims="-0.3;2.2";;
				# FEA4)	backcolor="blue"
						# trackcolor="blue"
						# fillcolor1="blue"
						# fillcolor2="blue"
						# lims="-2.5;5.5";;
						# # lims="-2.4;5.5";;
						# # lims="-1.2;5.1";;
				# TU1A)	backcolor="darkgreen"
						# trackcolor="darkgreen"
						# fillcolor1="darkgreen"
						# fillcolor2="darkgreen"
						# lims="-2.5;5.5";;
						# # lims="-2.4;5.5";;
						# # lims="-2.2;4.4";;
				# GT1)	backcolor="orange"
						# trackcolor="orange"
						# fillcolor1="orange"
						# fillcolor2="orange"
						# lims="-2.5;4";;
						# # lims="-1;3";;
				# HDZIV6)	backcolor="purple"
						# trackcolor="purple"
						# fillcolor1="purple"
						# fillcolor2="purple"
						# lims="-2.5;4";;
						# # lims="-2.4;5.5";;
						# # lims="-2.4;2.5";;
				# KN1)	backcolor="brown"
						# trackcolor="brown"
						# fillcolor1="brown"
						# fillcolor2="brown"
						# lims="-2.5;5.5";;
						# # lims="-2.4;5.5";;
						# # lims="-1.3;5.5";;
			# esac
			# bw="${pathtobw}/${line}_${name}_merged.bw"
			# printf "${line}\t${name}\t${bw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\tno\t${lims}\n" >> browser/data/${nameout}_samplefile.txt
		# fi
	# done < ${samplefile}
	
	# rm -f browser/data/${nameout}_loci.txt
	# while read short long GID
	# do
		# chr=$(grep "${GID}" manuscript/data/B73_v5_all_genes.bed | awk '{print $1}')
		# enstart=$(grep "${GID}" manuscript/data/B73_v5_all_genes.bed | awk '{print $2}')
		# enend=$(grep "${GID}" manuscript/data/B73_v5_all_genes.bed | awk '{print $3}')
		# enwidth=$((enend-enstart))
		# case "${short}" in
			# tb1)	bef="80000" #80000 ### 120,000
					# aft="20000";; #50000
			# fea4)	bef="60000" #10000 ### 30,000
					# aft="40000";; #20000
			# gt1)	bef="70000" #15000 ### 15,500
					# aft="30000";; #500
			# ra1)	bef="60000" #60000 ### 120,000
					# aft="40000";; #60000
			# tu1)	bef="20000" #20000
					# aft="20000";; #20000
			# hdziv6)	bef="6000" #6000
					# aft="100";; #100
			# kn1)	bef="5000" #5000
					# aft="120000";; #120000
			# tga1)	bef="40000"
					# aft="80000";;
		# esac
		# bs=10
		# start=$((enstart-bef))
		# end=$((enend+aft))
		# locus="${chr}:${start}:${end}"
		# printf "${locus}\t${long}\t${bs}\t${enstart}\t${enwidth}\n" >> browser/data/${nameout}_loci.txt
	# done < manuscript/data/tf2_genes_B73.txt
		
	# cd browser/
	# printf "\nRunning plot_browser_input_files script on ${nameout}\n"
	# qsub -sync y ../Plot_browser_input_files.sh -l B73 -r data/${nameout}_loci.txt -f data/${nameout}_samplefile.txt
	# cd ..
# done


####################################################################
## To do the upset plot for TSS clusters (Fig3c)
####################################################################

# line="B73"
# ref="B73_v5"
# analysisname="B73_on_B73_v5_all_genes"
# nameout="Fig2e_TSS"
# pathtss="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/RNA/TSS"
# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/B73_v5"
# regionfile="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/shRNA/tracks/B73_v5_all_genes.bed"

# i=0
# while read TEtype
# do
	# if [[ ${i} -eq 0 ]]; then
		# TEtypestring="${TEtype}"
	# else
		# TEtypestring="${TEtypestring},${TEtype}"
	# fi
	# i=$((i+1))
# done < combined/TSS/${ref}_TE_types.txt

# rm -f manuscript/data/tmp_${nameout}.bed
# for tissue in ears cn endosperm roots pollen
# do
	# printf "\nMaking tss file for ${tissue}\n"
	# awk -v OFS="\t" 'BEGIN {a=1} ($0 !~ /^#/) && ($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/) {print $1,$2,$3,"TSS_"a,"."; a+=1}' ${pathtss}/idr_${line}_${tissue}_RAMPAGE.narrowPeak | bedtools sort -g ${ref_dir}/chrom.sizes > manuscript/data/${line}_${tissue}_TSS.bed
	# awk -v OFS="\t" -v t=${tissue} '{print $1,$2,$3,t}' manuscript/data/${line}_${tissue}_TSS.bed | sort -k1,1 -k2,2n -u >> manuscript/data/tmp_${nameout}.bed
	# printf "\nGetting closest gene for ${tissue}\n"	
	# bedtools closest -a manuscript/data/${line}_${tissue}_TSS.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref -t first | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > manuscript/data/temp_${nameout}.bed
	# printf "\nGrouping based on distance for ${tissue}\n"
	# awk -v OFS="\t" '{if ($5<-2000) {d="Intergenic"} else if ($5<0) {d="Terminator"} else if ($5==0) {d="Gene_body"} else if ($5>2000) {d="Intergenic"} else {d="Promoter"} print $0,d}' manuscript/data/temp_${nameout}.bed > manuscript/data/temp2_${nameout}.bed
	# printf "\nIntersecting TE for ${tissue}\n"
	# bedtools intersect -a manuscript/data/temp2_${nameout}.bed -b combined/TSS/${ref}_all_tes.bed -loj | awk -v OFS="\t" -v t=${tissue} -v l=${line} '{if ($13==".") print l,t,l"_"t"_"$4,$9,"No",$9,$9; else if ($9 == "Intergenic") print l,t,l"_"t"_"$4,$9,$13,$13,$13; else print l,t,l"_"t"_"$4,$9,$13,$13,$13"_in_"$9}' > manuscript/data/temp_tss_in_genes_and_tes_${tissue}_${nameout}.bed
	# awk -v OFS="\t" 'BEGIN {a=""} {if ($3!=a) print; a=$3}' manuscript/data/temp_tss_in_genes_and_tes_${tissue}_${nameout}.bed > manuscript/data/tss_in_genes_and_tes_${tissue}_${nameout}.bed
	# rm -f manuscript/data/temp*_${nameout}.bed
# done
# printf "Line\tTissue\tTSS_ID\tGene\tTE\tLabel\tLabelcombined\n" > manuscript/data/Table_tissues_${nameout}.txt
# cat manuscript/data/tss_in_genes_and_tes_*_${nameout}.bed >> manuscript/data/Table_tissues_${nameout}.txt
# printf "\nPreparing merged tss file for ${nameout}\n"
# sort -k1,1 -k2,2n manuscript/data/tmp_${nameout}.bed > manuscript/data/tmp2_${nameout}.bed
# bedtools merge -i manuscript/data/tmp2_${nameout}.bed -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{print $1,$2,$3,"TSS_"NR,$4}'> manuscript/data/tmp3_${nameout}.bed
# printf "\nGetting closest gene for tss in ${nameout}\n"
# bedtools closest -a manuscript/data/tmp3_${nameout}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > manuscript/data/tmp4_${nameout}.bed
# printf "\nGrouping based on distance\n"
# awk -v OFS="\t" '{if ($5<-2000) {d="Intergenic"} else if ($5<0) {d="Terminator"} else if ($5==0) {d="Gene_body"} else if ($5>2000) {d="Intergenic"} else {d="Promoter"} print $0,d}' manuscript/data/tmp4_${nameout}.bed > manuscript/data/tmp5_${nameout}.bed
# printf "\nIntersecting with TEs\n"
# bedtools intersect -a manuscript/data/tmp5_${nameout}.bed -b combined/TSS/${ref}_all_tes.bed -loj | awk -v OFS="\t" -v l=${line} 'BEGIN {printf "Line\tTSS_ID\tGene\tTE\tLabel\tLabelcombined\tGID\tTissues\tPeak_coordinates\n"} {if ($13==".") print l,$4,$9,"No",$9,$9,$8,$7,$1":"$2"-"$3; else if ($9 == "Intergenic") print l,$4,$9,$13,$13,$13,$8,$7,$1":"$2"-"$3; else print l,$4,$9,$13,$13,$13"_in_"$9,$8,$7,$1":"$2"-"$3}' > manuscript/data/tmp_all_tss_in_genes_and_tes_shRNA_${nameout}.bed
# awk -v OFS="\t" 'BEGIN {a=""} {if ($2!=a) print; a=$2}' manuscript/data/tmp_all_tss_in_genes_and_tes_shRNA_${nameout}.bed > manuscript/data/all_tss_in_genes_and_tes_${nameout}.bed
# rm -f manuscript/data/tmp*_${nameout}.bed
# #### To create a matrix of peak presence in each sample
# printf "\nCreating matrix file for ${nameout}\n"
# for tissue in ears cn endosperm roots pollen
# do
	# printf "${tissue}\n" > manuscript/data/temp_col_${nameout}_${tissue}.txt
	# awk -v OFS="\t" -v t=${tissue} 'NR>1 {if ($8 ~ t) print "1"; else print "0"}' manuscript/data/all_tss_in_genes_and_tes_${nameout}.bed >> manuscript/data/temp_col_${nameout}_${tissue}.txt
# done
# awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7}' manuscript/data/all_tss_in_genes_and_tes_${nameout}.bed > manuscript/data/temp_col_${nameout}_AAA.txt
# paste manuscript/data/temp_col_${nameout}_*.txt | uniq > manuscript/data/matrix_upset_${nameout}.txt
# rm -f manuscript/data/temp_col_${nameout}_*.txt
# ### To make an Upset plot highlighting peaks in gene bodies
# printf "\nCreating Distirbution and Upset plot for shRNA clusters in ${nameout} with R version:\n"
# R --version
# Rscript --vanilla MaizeCode_R_TSS_distribution_upset.r ${nameout} ${TEtypestring} manuscript/data/Table_tissues_${nameout}.txt manuscript/data/matrix_upset_${nameout}.txt

####################################################################
## To do the upset plot for shRNA clusters (Fig3d)
####################################################################

# line="B73"
# ref="B73_v5"
# analysisname="B73_on_B73_v5_all_genes"
# nameout="Fig2f_shRNA_cluster"
# pathshrna="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/shRNA/mapped"
# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/B73_v5"
# regionfile="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/shRNA/tracks/B73_v5_all_genes.bed"

# i=0
# while read TEtype
# do
	# if [[ ${i} -eq 0 ]]; then
		# TEtypestring="${TEtype}"
	# else
		# TEtypestring="${TEtypestring},${TEtype}"
	# fi
	# i=$((i+1))
# done < combined/TSS/${ref}_TE_types.txt

# rm -f manuscript/data/tmp_shRNA_clusters_${nameout}.bed
# for tissue in ears cn endosperm roots pollen
# do
	# printf "\nMaking shRNA cluster file for ${tissue}\n"
	# awk -v OFS="\t" 'BEGIN {a=1} ($0 !~ /^#/) && ($1~/^[0-9]/ || $1~/^chr[0-9]/ || $1~/^Chr[0-9]/) {print $1,$4,$5,$3"_"a,"."; a+=1}' ${pathshrna}/${line}_${tissue}_shRNA/ShortStack_All.gff3 | bedtools sort -g ${ref_dir}/chrom.sizes > manuscript/data/${line}_${tissue}_shRNA_clusters.bed
	# awk -v OFS="\t" -v t=${tissue} '{print $1,$2,$3,t}' manuscript/data/${line}_${tissue}_shRNA_clusters.bed | sort -k1,1 -k2,2n -u >> manuscript/data/tmp_shRNA_clusters_${nameout}.bed
	# printf "\nGetting closest gene for ${tissue}\n"	
	# bedtools closest -a manuscript/data/${line}_${tissue}_shRNA_clusters.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref -t first | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > manuscript/data/temp_clusters_${nameout}.bed
	# printf "\nGrouping based on distance for ${tissue}\n"
	# awk -v OFS="\t" '{if ($5<-2000) {d="Intergenic"} else if ($5<0) {d="Terminator"} else if ($5==0) {d="Gene_body"} else if ($5>2000) {d="Intergenic"} else {d="Promoter"} print $0,d}' manuscript/data/temp_clusters_${nameout}.bed > manuscript/data/temp2_clusters_${nameout}.bed
	# printf "\nIntersecting TE for ${tissue}\n"
	# bedtools intersect -a manuscript/data/temp2_clusters_${nameout}.bed -b combined/TSS/${ref}_all_tes.bed -loj | awk -v OFS="\t" -v t=${tissue} -v l=${line} '{if ($13==".") print l,t,l"_"t"_"$4,$9,"No",$9,$9; else if ($9 == "Intergenic") print l,t,l"_"t"_"$4,$9,$13,$13,$13; else print l,t,l"_"t"_"$4,$9,$13,$13,$13"_in_"$9}' > manuscript/data/temp_clusters_in_genes_and_tes_${tissue}_${nameout}.bed
	# awk -v OFS="\t" 'BEGIN {a=""} {if ($3!=a) print; a=$3}' manuscript/data/temp_clusters_in_genes_and_tes_${tissue}_${nameout}.bed > manuscript/data/clusters_in_genes_and_tes_${tissue}_${nameout}.bed
	# rm -f manuscript/data/temp*_clusters_${nameout}.bed
# done
# printf "Line\tTissue\tCluster_ID\tGene\tTE\tLabel\tLabelcombined\n" > manuscript/data/Table_shRNA_clusters_tissues_${nameout}.txt
# cat manuscript/data/clusters_in_genes_and_tes_*_${nameout}.bed >> manuscript/data/Table_shRNA_clusters_tissues_${nameout}.txt
# printf "\nPreparing merged shRNA cluster file for ${nameout}\n"
# sort -k1,1 -k2,2n manuscript/data/tmp_shRNA_clusters_${nameout}.bed > manuscript/data/tmp2_shRNA_clusters_${nameout}.bed
# bedtools merge -i manuscript/data/tmp2_shRNA_clusters_${nameout}.bed -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{print $1,$2,$3,"Cluster_"NR,$4}'> manuscript/data/tmp3_shRNA_clusters_${nameout}.bed
# printf "\nGetting closest gene for clusters in ${nameout}\n"
# bedtools closest -a manuscript/data/tmp3_shRNA_clusters_${nameout}.bed -b ${regionfile} -g ${ref_dir}/chrom.sizes -D ref | awk -v OFS="\t" '{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}' | awk -F"[:=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$9}' > manuscript/data/tmp4_shRNA_clusters_${nameout}.bed
# printf "\nGrouping based on distance\n"
# awk -v OFS="\t" '{if ($5<-2000) {d="Intergenic"} else if ($5<0) {d="Terminator"} else if ($5==0) {d="Gene_body"} else if ($5>2000) {d="Intergenic"} else {d="Promoter"} print $0,d}' manuscript/data/tmp4_shRNA_clusters_${nameout}.bed > manuscript/data/tmp5_shRNA_clusters_${nameout}.bed
# printf "\nIntersecting with TEs\n"
# bedtools intersect -a manuscript/data/tmp5_shRNA_clusters_${nameout}.bed -b combined/TSS/${ref}_all_tes.bed -loj | awk -v OFS="\t" -v l=${line} 'BEGIN {printf "Line\tCluster_ID\tGene\tTE\tLabel\tLabelcombined\tGID\tTissues\tPeak_coordinates\n"} {if ($13==".") print l,$4,$9,"No",$9,$9,$8,$7,$1":"$2"-"$3; else if ($9 == "Intergenic") print l,$4,$9,$13,$13,$13,$8,$7,$1":"$2"-"$3; else print l,$4,$9,$13,$13,$13"_in_"$9,$8,$7,$1":"$2"-"$3}' > manuscript/data/temp_all_clusters_in_genes_and_tes_shRNA_${nameout}.bed
# awk -v OFS="\t" 'BEGIN {a=""} {if ($2!=a) print; a=$2}' manuscript/data/temp_all_clusters_in_genes_and_tes_shRNA_${nameout}.bed > manuscript/data/all_shRNA_clusters_in_genes_and_tes_${nameout}.bed
# rm -f manuscript/data/tmp*_shRNA_clusters_${nameout}.bed
# #### To create a matrix of peak presence in each sample
# printf "\nCreating matrix file for ${nameout}\n"
# for tissue in ears cn endosperm roots pollen
# do
	# printf "${tissue}\n" > manuscript/data/temp_col_clusters_${nameout}_${tissue}.txt
	# awk -v OFS="\t" -v t=${tissue} 'NR>1 {if ($8 ~ t) print "1"; else print "0"}' manuscript/data/all_shRNA_clusters_in_genes_and_tes_${nameout}.bed >> manuscript/data/temp_col_clusters_${nameout}_${tissue}.txt
# done
# awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7}' manuscript/data/all_shRNA_clusters_in_genes_and_tes_${nameout}.bed > manuscript/data/temp_col_clusters_${nameout}_AAA.txt
# paste manuscript/data/temp_col_clusters_${nameout}_*.txt | uniq > manuscript/data/matrix_upset_shRNA_clusters_${nameout}.txt
# rm -f manuscript/data/temp_col_clusters_${nameout}_*.txt
# ### To make an Upset plot highlighting peaks in gene bodies
# printf "\nCreating Distirbution and Upset plot for shRNA clusters in ${nameout} with R version:\n"
# R --version
# Rscript --vanilla MaizeCode_R_shRNA_distribution_upset.r ${nameout} ${TEtypestring} manuscript/data/Table_shRNA_clusters_tissues_${nameout}.txt manuscript/data/matrix_upset_shRNA_clusters_${nameout}.txt

############################################################################################################################
# # #### To get the different groups of enhancers (with and without H3K4me1, with and without expression) and plot the heatmaps and profiles (Fig 4a,c,d,e,f and Sup Fig 6a,b)
############################################################################################################################

# analysisname="B73_on_B73_v5_all_genes"
# pathtobwchip="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/ChIP/tracks"
# pathtobwrna="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/RNA/tracks"
# pathtobwshrna="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/shRNA/tracks"
# pathtobamshrna="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/shRNA/mapped"
# pathtobwmnase="/grid/martienssen/home/jcahn/norepl/projects/MNase/ChIP/tracks"
# pathtofiles="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode"

# line="B73"
# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/B73_v5"
# ref="B73_v5"
# colorsmap=("Reds" "Blues" "YlOrBr" "Purples" "Greys" "Reds" "Blues" "Greys" "Greys")
# # colorsprofile=("#5546EE" "#33A9ED" "#33CBED" "#ED4433" "#ED8A33" "#EDB533")
# colorsprofile=("#A9A9A9" "#5546EE" "#ED4433" "#ED8A33" "#EDB533")

# # grep -w "knob" combined/TSS/B73_v5_all_tes.bed > manuscript/data/blacklist_regions.bed
# # grep "centromeric" combined/TSS/B73_v5_all_tes.bed >> manuscript/data/blacklist_regions.bed

# # path="/grid/martienssen/home/jcahn/nlsas/projects/methyl_maize"
# # mcbw=()
# # mcsamplelab=()
# # mcmins=()
# # mcmaxs=()
# # for context in CG CHG CHH
# # do
	# # case "${context}" in
		# # CG)		maxi1=100;;
		# # CHG)	maxi1=80;;
		# # CHH)	maxi1=8;;
	# # esac
	# # for tis in seedlings
	# # do
		# # for rep in Rep1
		# # do
			# # name="B73_${tis}_mC_${rep}_${context}"
			# # mcbw+=("${path}/mC/tracks/${name}.bw")
			# # mcsamplelab+=("${tis}_${context}")
			# # mcmins+=("0")
			# # mcmaxs+=("${maxi1}")
		# # done
	# # done
# # done

# # for tissue in cn ears roots endosperm
# for tissue in cn
# do
	# # # To prepare the RNA, RAMPAGE, shRNA and siRNA tracks, if needed
	# # for rna in RNAseq RAMPAGE
	# # do
		# # if [ ! -e ${pathtobwrna}/${line}_${tissue}_${rna}_merged_sum.bw ]; then
			# # printf "merging strands for ${rna}\n"
			# # bigWigMerge ${pathtobwrna}/${line}_${tissue}_${rna}_merged_plus.bw ${pathtobwrna}/${line}_${tissue}_${rna}_merged_minus.bw ${pathtobwrna}/${line}_${tissue}_${rna}_merged_sum.bg
			# # sort -k1,1 -k2,2n ${pathtobwrna}/${line}_${tissue}_${rna}_merged_sum.bg > ${pathtobwrna}/${line}_${tissue}_${rna}_merged_sum_sorted.bg
			# # bedGraphToBigWig ${pathtobwrna}/${line}_${tissue}_${rna}_merged_sum_sorted.bg ${ref_dir}/chrom.sizes ${pathtobwrna}/${line}_${tissue}_${rna}_merged_sum.bw
			# # rm -f ${pathtobwrna}/${line}_${tissue}_${rna}*.bg
		# # fi
	# # done
	
	# # if [ ! -e ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_sum.bw ]; then
		# # printf "Filtering ${line}_${tissue}_Rep1 shRNA files for reads >30nt\n"
		# # samtools view -@ ${threads} -h ${pathtobamshrna}/${line}_${tissue}_shRNA_Rep1/filtered_${line}_${tissue}_shRNA_Rep1.bam | awk 'length($10) >= 30 || $1 ~ /^@/' | samtools view -bS - > ${pathtobamshrna}/temp_${line}_${tissue}_shortRNA_Rep1.bam
		# # printf "Filtering ${line}_${tissue}_Rep2 shRNA files for reads >30nt\n"
		# # samtools view -@ ${threads} -h ${pathtobamshrna}/${line}_${tissue}_shRNA_Rep2/filtered_${line}_${tissue}_shRNA_Rep2.bam | awk 'length($10) >= 30 || $1 ~ /^@/' | samtools view -bS - > ${pathtobamshrna}/temp_${line}_${tissue}_shortRNA_Rep2.bam
		# # printf "Merging ${line}_${tissue} Replicates\n"
		# # samtools merge -@ ${threads} ${pathtobamshrna}/temp_${line}_${tissue}_shortRNA_merged.bam ${pathtobamshrna}/temp_${line}_${tissue}_shortRNA_Rep*.bam
		# # samtools sort -@ ${threads} -o ${pathtobamshrna}/${line}_${tissue}_shortRNA_merged.bam ${pathtobamshrna}/temp_${line}_${tissue}_shortRNA_merged.bam
		# # rm -f ${pathtobamshrna}/temp_${line}_${tissue}_shortRNA*.bam
		# # samtools index -@ ${threads} ${pathtobamshrna}/${line}_${tissue}_shortRNA_merged.bam
		# # printf "Getting stranded coverage for ${line}_${tissue}\n"
		# # bamCoverage --filterRNAstrand forward -bs 1 -p ${threads} --normalizeUsing CPM -b ${pathtobamshrna}/${line}_${tissue}_shortRNA_merged.bam -o ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_plus.bw
		# # bamCoverage --filterRNAstrand reverse -bs 1 -p ${threads} --normalizeUsing CPM -b ${pathtobamshrna}/${line}_${tissue}_shortRNA_merged.bam -o ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_minus.bw
		# # printf "Merging strands for ${line}_${tissue} shortRNA\n"
		# # bigWigMerge ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_plus.bw ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_minus.bw ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_sum.bg
		# # sort -k1,1 -k2,2n ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_sum.bg > ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_sum_sorted.bg
		# # bedGraphToBigWig ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_sum_sorted.bg ${ref_dir}/chrom.sizes ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_sum.bw
		# # rm -f ${pathtobwshrna}/${line}_${tissue}_shortRNA*.bg
	# # fi
	
	# # if [ ! -e ${pathtobwshrna}/${line}_${tissue}_siRNA_merged_sum.bw ]; then
		# # printf "Merging strands for ${line}_${tissue} siRNA\n"
		# # bigWigMerge ${pathtobwshrna}/${line}_${tissue}_shRNA_merged_plus.bw ${pathtobwshrna}/${line}_${tissue}_shRNA_merged_minus.bw ${pathtobwshrna}/${line}_${tissue}_siRNA_merged_sum.bg
		# # sort -k1,1 -k2,2n ${pathtobwshrna}/${line}_${tissue}_siRNA_merged_sum.bg > ${pathtobwshrna}/${line}_${tissue}_siRNA_merged_sum_sorted.bg
		# # bedGraphToBigWig ${pathtobwshrna}/${line}_${tissue}_siRNA_merged_sum_sorted.bg ${ref_dir}/chrom.sizes ${pathtobwshrna}/${line}_${tissue}_siRNA_merged_sum.bw
		# # rm -f ${pathtobwshrna}/${line}_${tissue}_siRNA*.bg
	# # fi
	
	# # if [ ! -e ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum.bw ]; then
		# # printf "Filtering ${line}_${tissue}_Rep1 shRNA files for reads =24nt\n"
		# # samtools view -@ ${threads} -h ${pathtobamshrna}/${line}_${tissue}_shRNA_Rep1/filtered_${line}_${tissue}_shRNA_Rep1.bam | awk 'length($10) == 24 || $1 ~ /^@/' | samtools view -bS - > ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_Rep1.bam
		# # printf "Filtering ${line}_${tissue}_Rep2 shRNA files for reads =24nt\n"
		# # samtools view -@ ${threads} -h ${pathtobamshrna}/${line}_${tissue}_shRNA_Rep2/filtered_${line}_${tissue}_shRNA_Rep2.bam | awk 'length($10) == 24 || $1 ~ /^@/' | samtools view -bS - > ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_Rep2.bam
		# # printf "Merging ${line}_${tissue} Replicates\n"
		# # samtools merge -@ ${threads} ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_merged.bam ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_Rep*.bam
		# # samtools sort -@ ${threads} -o ${pathtobamshrna}/${line}_${tissue}_24RNA_merged.bam ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_merged.bam
		# # rm -f ${pathtobamshrna}/temp_${line}_${tissue}_24RNA*.bam
		# # samtools index -@ ${threads} ${pathtobamshrna}/${line}_${tissue}_24RNA_merged.bam
		# # printf "Getting stranded coverage for ${line}_${tissue}\n"
		# # bamCoverage --filterRNAstrand forward -bs 1 -p ${threads} --normalizeUsing CPM -b ${pathtobamshrna}/${line}_${tissue}_24RNA_merged.bam -o ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_plus.bw
		# # bamCoverage --filterRNAstrand reverse -bs 1 -p ${threads} --normalizeUsing CPM -b ${pathtobamshrna}/${line}_${tissue}_24RNA_merged.bam -o ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_minus.bw
		# # printf "Merging strands for ${line}_${tissue} 24RNA\n"
		# # bigWigMerge ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_plus.bw ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_minus.bw ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum.bg
		# # sort -k1,1 -k2,2n ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum.bg > ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum_sorted.bg
		# # bedGraphToBigWig ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum_sorted.bg ${ref_dir}/chrom.sizes ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum.bw
		# # rm -f ${pathtobwshrna}/${line}_${tissue}_24RNA*.bg
	# # fi

	# # printf "getting regions for ${nameout}\n"
	# # awk -v OFS="\t" 'NR>1' ${pathtofiles}/combined/peaks/complete_enhancers_distal_upstream_B73_${tissue}_${analysisname}.txt > manuscript/data/temp_enhancers_distal_B73_${tissue}.txt
	# # awk -v OFS="\t" 'NR>1' ${pathtofiles}/combined/peaks/complete_enhancers_distal_downstream_B73_${tissue}_${analysisname}.txt >> manuscript/data/temp_enhancers_distal_B73_${tissue}.txt
	# # bedtools sort -g ${ref_dir}/chrom.sizes -i manuscript/data/temp_enhancers_distal_B73_${tissue}.txt > manuscript/data/temp2_enhancers_distal_B73_${tissue}.txt
	# # bedtools intersect -v -wa -a manuscript/data/temp2_enhancers_distal_B73_${tissue}.txt -b manuscript/data/blacklist_regions.bed > manuscript/data/enhancers_distal_B73_${tissue}.txt
	# # rm -f manuscript/data/temp*_enhancers_distal_B73_${tissue}.txt
	# # bedtools closest -a manuscript/data/enhancers_distal_B73_${tissue}.txt -b ${pathtofiles}/ChIP/peaks/best_peaks_B73_${tissue}_H3K4me1.bed -D ref -t first -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" -v t=${tissue} '{if ($26>1000 || $26<-1000 ) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14 > "manuscript/data/temp_results_"t"_without_H3K4me1.txt"; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14 > "manuscript/data/temp_results_"t"_with_H3K4me1.txt"}'
		
	# # awk -v OFS="\t" -v t=${tissue} '{if ($9>0 && $10>0) print $4 > "manuscript/data/temp_"t"_with_k4_expdouble.txt"; else if ($9>0 || $10>0) print $4 > "manuscript/data/temp_"t"_with_k4_expsingle.txt"; else print $4 > "manuscript/data/temp_"t"_with_k4_noexp.txt";}' manuscript/data/temp_results_${tissue}_with_H3K4me1.txt
	# # awk -v OFS="\t" -v t=${tissue} '{if ($9>0 && $10>0) print $4 > "manuscript/data/temp_"t"_without_k4_expdouble.txt"; else if ($9>0 || $10>0) print $4 > "manuscript/data/temp_"t"_without_k4_expsingle.txt"; else print $4 > "manuscript/data/temp_"t"_without_k4_noexp.txt";}' manuscript/data/temp_results_${tissue}_without_H3K4me1.txt
	
	# # nameout="Fig4_B73_${tissue}_distal_enhancers_all"
	# # nameout="Fig4_B73_${tissue}_distal_enhancers_without"

	# labels="${tissue}_H3K27ac ${tissue}_H3K4me1 ${tissue}_H3K4me3 ${tissue}_RNAseq ${tissue}_RAMPAGE ${tissue}_shortRNA ${tissue}_siRNA ${tissue}_24RNA ${tissue}_MNase conservation_score"
	# # labels="${tissue}_H3K27ac ${tissue}_H3K4me1 ${tissue}_H3K4me3 ${tissue}_RNAseq ${tissue}_RAMPAGE ${tissue}_shortRNA ${tissue}_siRNA"
	
	# region_label=()
	# region_list_ref=()
	# region_list_whole=()
	
	# # cat manuscript/data/B73_v5_${tissue}_with_k4_expdouble_summits.bed manuscript/data/B73_v5_${tissue}_with_k4_expsingle_summits.bed manuscript/data/B73_v5_${tissue}_with_k4_noexp_summits.bed | sort -k5,5nr > manuscript/data/B73_v5_${tissue}_with_k4_summits.bed
	# # cat manuscript/data/B73_v5_${tissue}_with_k4_expdouble_enhancers.bed manuscript/data/B73_v5_${tissue}_with_k4_expsingle_enhancers.bed manuscript/data/B73_v5_${tissue}_with_k4_noexp_enhancers.bed | sort -k5,5nr > manuscript/data/B73_v5_${tissue}_with_k4_enhancers.bed
	# # for type in with_k4_expdouble with_k4_expsingle with_k4_noexp without_k4_expdouble without_k4_expsingle without_k4_noexp
	# # for type in with_k4 without_k4_expdouble without_k4_expsingle without_k4_noexp
	# for type in shuffled with_k4 without_k4_expdouble without_k4_expsingle without_k4_noexp
	# do
		# # printf "getting ${type} enhancers for ${tissue}\n"
		# # grep -w -f manuscript/data/temp_${tissue}_${type}.txt ${pathtofiles}/ChIP/peaks/best_peaks_B73_${tissue}_H3K27ac.bed | awk -v OFS="\t" '{print $1,$2+$10,$2+$10+1,$4,$5,$6}' > manuscript/data/temp2_${tissue}_${type}.txt
		# # sort -k5,5nr manuscript/data/temp2_${tissue}_${type}.txt > manuscript/data/B73_v5_${tissue}_${type}_summits.bed
		# # grep -w -f manuscript/data/temp_${tissue}_${type}.txt ${pathtofiles}/ChIP/peaks/best_peaks_B73_${tissue}_H3K27ac.bed | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' > manuscript/data/temp3_${tissue}_${type}.txt
		# # sort -k5,5nr manuscript/data/temp3_${tissue}_${type}.txt > manuscript/data/B73_v5_${tissue}_${type}_enhancers.bed

		# # rm -f manuscript/data/temp*_${tissue}_${type}.txt
		# n=$(wc -l manuscript/data/B73_v5_${tissue}_${type}_summits.bed | awk '{print $1}')
		# region_label+=("${type}($n)")
		# region_list_ref+=("manuscript/data/B73_v5_${tissue}_${type}_summits.bed")
		# region_list_whole+=("manuscript/data/B73_v5_${tissue}_${type}_enhancers.bed")
	# done
	
	# bw_list="${pathtobwchip}/${line}_${tissue}_H3K27ac_merged.bw ${pathtobwchip}/${line}_${tissue}_H3K4me1_merged.bw ${pathtobwchip}/${line}_${tissue}_H3K4me3_merged.bw ${pathtobwrna}/${line}_${tissue}_RNAseq_merged_sum.bw ${pathtobwrna}/${line}_${tissue}_RAMPAGE_merged_sum.bw ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_sum.bw ${pathtobwshrna}/${line}_${tissue}_siRNA_merged_sum.bw ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum.bw ${pathtobwmnase}/${line}_${tissue}_MNase_merged.bw Armin_conservation_pan_andropogoneae/B73_v5_conservation_scores.bw ${pathtofiles}/combined/tracks/B73_v5_all_genes.bw ${pathtofiles}/combined/tracks/B73_v5_all_tes.bw"
	# # bw_list="${pathtobwchip}/${line}_${tissue}_H3K27ac_merged.bw ${pathtobwchip}/${line}_${tissue}_H3K4me1_merged.bw ${pathtobwchip}/${line}_${tissue}_H3K4me3_merged.bw ${pathtobwrna}/${line}_${tissue}_RNAseq_merged_sum.bw ${pathtobwrna}/${line}_${tissue}_RAMPAGE_merged_sum.bw ${pathtobwshrna}/${line}_${tissue}_shortRNA_merged_sum.bw ${pathtobwshrna}/${line}_${tissue}_siRNA_merged_sum.bw ${pathtofiles}/combined/tracks/B73_v5_all_genes.bw ${pathtofiles}/combined/tracks/B73_v5_all_tes.bw"
	
	# #### To make heatmap on summits
	
	# # printf "\nComputing scale-regions matrix for ${nameout}\n"
	# # computeMatrix reference-point --referencePoint "TSS" --missingDataAsZero --skipZeros -R ${region_list_ref[@]} -S ${bw_list} -bs 50 -b 5000 -a 5000 -p ${threads} -o manuscript/matrix/${nameout}.gz --quiet

	# # printf "\nGetting scales for ${nameout}\n"
	# # computeMatrixOperations dataRange -m manuscript/matrix/${nameout}.gz > manuscript/matrix/values_z_${nameout}.txt
	# # mins=()
	# # maxs=()
	# # for sample in ${labels[@]}
	# # do
		# # mini=$(grep ${sample} manuscript/matrix/values_z_${nameout}.txt | awk -v OFMT="%f" '{print $5}')
		# # maxi=$(grep ${sample} manuscript/matrix/values_z_${nameout}.txt | awk -v OFMT="%f" '{print $6}')
		# # test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		# # if [[ ${test} == "yes" ]]; then
			# # mins+=("0")
			# # maxs+=("0.00001")
		# # else
			# # mins+=("${mini}")
			# # maxs+=("${maxi}")
		# # fi
	# # done
	# # mins+=("0" "0")
	# # maxs+=("1" "1")

	# # plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_profile.pdf --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --averageType mean --outFileNameData manuscript/matrix/values_y_${nameout}.txt
	# # ymins=()
	# # ymaxs=()
	# # for sample in ${labels[@]}
	# # do
		# # ymini=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk -v OFMT="%f" 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
		# # ymaxi=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk -v OFMT="%f" 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
		# # test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		# # if [[ ${test} == "yes" ]]; then
			# # ymins+=("0")
			# # ymaxs+=("0.00001")
		# # else
			# # ymins+=("${ymini}")
			# # ymaxs+=("${ymaxi}")
		# # fi
	# # done
	# # ymins+=("0" "0")
	# # ymaxs+=("1" "1")
	# # printf "\nPlotting heatmap for ${nameout} scaling by sample\n"
	# # plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_ref.pdf --sortRegions keep --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --refPointLabel "enhancer" --colorMap ${colorsmap[@]} --outFileSortedRegions manuscript/data/${tissue}_sorted_ref_type_enhancers.txt
	# # printf "\nPlotting heatmap for ${nameout} scaling by sample\n"
	# # plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_ref_heatmap.pdf --sortRegions keep --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --refPointLabel "enhancer" --whatToShow "heatmap only" --colorMap ${colorsmap[@]}
	# # plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_ref_heatmap_bar.pdf --sortRegions keep --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --refPointLabel "enhancer" --whatToShow "heatmap and colorbar" --colorMap ${colorsmap[@]}
	# # printf "\nPlotting mean profile for ${nameout} scaling by sample\n"
	# # plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_ref_profile_mean.pdf --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --averageType mean --yMin ${ymins[@]} --yMax ${ymaxs[@]} --refPointLabel "enhancer" --colors ${colorsprofile[@]}
	
	# # plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_ref_profile.pdf --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --averageType median --outFileNameData manuscript/matrix/values_y_${nameout}.txt
	# # ymins=()
	# # ymaxs=()
	# # for sample in ${labels[@]}
	# # do
	 	# # ymini=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk -v OFMT="%f" 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
		# # ymaxi=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk -v OFMT="%f" 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
		# # test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		# # if [[ ${test} == "yes" ]]; then
			# # ymins+=("0")
			# # ymaxs+=("0.00001")
		# # else
			# # ymins+=("${ymini}")
			# # ymaxs+=("${ymaxi}")
		# # fi
	# # done
	# # ymins+=("0" "0")
	# # ymaxs+=("1" "1")
	# # printf "\nPlotting median profile for ${nameout} scaling by sample\n"
	# # plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_ref_profile_median.pdf --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --averageType median --yMin ${ymins[@]} --yMax ${ymaxs[@]} --refPointLabel "enhancer" --colors ${colorsprofile[@]}
	
	# # rm -f manuscript/data/temp*
	# # rm -f manuscript/matrix/temp*
	# # rm -f manuscript/plots/${nameout}_ref_profile.pdf
	# # rm -f manuscript/data/distal_${tissue}_*.txt
	
	# # printf "Computing mC matrix ${tissue}\n"
	# # computeMatrix reference-point -R manuscript/data/${tissue}_sorted_ref_type_enhancers.txt -S ${mcbw[@]} -bs 50 -b 5000 -a 5000 --sortRegions keep -p ${threads} -o manuscript/matrix/matrix_mC_${tissue}.gz
	# # printf "Plotting heatmap ${tissue}\n"
	# # plotHeatmap -m manuscript/matrix/matrix_mC_${tissue}.gz -out manuscript/plots/${nameout}_ref_mC.pdf --sortRegions keep --samplesLabel ${mcsamplelab[@]} --regionsLabel ${region_label[@]} --colorMap 'Oranges' --missingDataColor 'grey' --interpolationMethod 'nearest' --yMin ${mcmins[@]} --yMax ${mcmaxs[@]} --zMin ${mcmins[@]} --zMax ${mcmaxs[@]} --refPointLabel "enhancer"
	# # plotHeatmap -m manuscript/matrix/matrix_mC_${tissue}.gz -out manuscript/plots/${nameout}_ref_heatmap_mC.pdf --sortRegions keep --samplesLabel ${mcsamplelab[@]} --regionsLabel ${region_label[@]} --colorMap 'Oranges' --missingDataColor 'grey' --interpolationMethod 'nearest' --yMin ${mcmins[@]} --yMax ${mcmaxs[@]} --zMin ${mcmins[@]} --zMax ${mcmaxs[@]} --refPointLabel "enhancer" --whatToShow "heatmap only"
	# # plotProfile -m manuscript/matrix/matrix_mC_${tissue}.gz -out manuscript/plots/${nameout}_ref_profile_mC_mean.pdf --samplesLabel ${mcsamplelab[@]} --regionsLabel ${region_label[@]} --numPlotsPerRow 3 --yMin ${mcmins[@]} --yMax ${mcmaxs[@]} --averageType mean --refPointLabel "enhancer" --colors ${colorsprofile[@]}
	# # plotProfile -m manuscript/matrix/matrix_mC_${tissue}.gz -out manuscript/plots/${nameout}_ref_profile_mC_median.pdf --samplesLabel ${mcsamplelab[@]} --regionsLabel ${region_label[@]} --numPlotsPerRow 3 --yMin ${mcmins[@]} --yMax ${mcmaxs[@]} --averageType median --refPointLabel "enhancer" --colors ${colorsprofile[@]}

	# #### To make heatmap on whole enhancers
	
	# printf "\nComputing scale-regions matrix for ${nameout}\n"
	# computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${region_list_whole[@]} -S ${bw_list} -bs 50 -b 5000 -a 5000 -m 1000 -p ${threads} -o manuscript/matrix/${nameout}.gz --quiet

	# # printf "\nGetting scales for ${nameout}\n"
	# # computeMatrixOperations dataRange -m manuscript/matrix/${nameout}.gz > manuscript/matrix/values_z_${nameout}.txt
	# # mins=()
	# # maxs=()
	# # for sample in ${labels[@]}
	# # do
		# # mini=$(grep ${sample} manuscript/matrix/values_z_${nameout}.txt | awk -v OFMT="%f" '{print $5}')
		# # maxi=$(grep ${sample} manuscript/matrix/values_z_${nameout}.txt | awk -v OFMT="%f" '{print $6}')
		# # test=$(awk -v a=${mini} -v b=${maxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		# # if [[ ${test} == "yes" ]]; then
			# # mins+=("0")
			# # maxs+=("0.00001")
		# # else
			# # mins+=("${mini}")
			# # maxs+=("${maxi}")
		# # fi
	# # done
	# # mins+=("0" "0")
	# # maxs+=("1" "1")

	# plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_profile.pdf --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --averageType mean --outFileNameData manuscript/matrix/values_y_${nameout}.txt
	# ymins=()
	# ymaxs=()
	# for sample in ${labels[@]}
	# do
		# ymini=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk -v OFMT="%f" 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
		# ymaxi=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk -v OFMT="%f" 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
		# test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		# if [[ ${test} == "yes" ]]; then
			# ymins+=("0")
			# ymaxs+=("0.00001")
		# else
			# ymins+=("${ymini}")
			# ymaxs+=("${ymaxi}")
		# fi
	# done
	# ymins+=("0" "0")
	# ymaxs+=("1" "1")
	# # printf "\nPlotting heatmap for ${nameout} scaling by sample\n"
	# # plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_scaled.pdf --sortRegions keep --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --xAxisLabel "enhancer" --colorMap ${colorsmap[@]} --outFileSortedRegions manuscript/data/${tissue}_sorted_scaled_type_enhancers.txt
	# # printf "\nPlotting heatmap for ${nameout} scaling by sample\n"
	# # plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_scaled_heatmap.pdf --sortRegions keep --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --xAxisLabel "enhancer" --whatToShow "heatmap only" --colorMap ${colorsmap[@]}
	# # plotHeatmap -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_scaled_heatmap_bar.pdf --sortRegions keep --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --zMin ${mins[@]} --zMax ${maxs[@]} --yMin ${ymins[@]} --yMax ${ymaxs[@]} --interpolationMethod 'bilinear' --xAxisLabel "enhancer" --whatToShow "heatmap and colorbar" --colorMap ${colorsmap[@]}
	# printf "\nPlotting mean profile for ${nameout} scaling by sample\n"
	# plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_scaled_profile_mean.pdf --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --averageType mean --yMin ${ymins[@]} --yMax ${ymaxs[@]} --colors ${colorsprofile[@]}
	
	# # plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_scaled_profile.pdf --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --averageType median --outFileNameData manuscript/matrix/values_y_${nameout}.txt
	# # ymins=()
	# # ymaxs=()
	# # for sample in ${labels[@]}
	# # do
	 	# # ymini=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk -v OFMT="%f" 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.8; print a}')
		# # ymaxi=$(grep ${sample} manuscript/matrix/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk -v OFMT="%f" 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
		# # test=$(awk -v a=${ymini} -v b=${ymaxi} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		# # if [[ ${test} == "yes" ]]; then
			# # ymins+=("0")
			# # ymaxs+=("0.00001")
		# # else
			# # ymins+=("${ymini}")
			# # ymaxs+=("${ymaxi}")
		# # fi
	# # done
	# # ymins+=("0" "0")
	# # ymaxs+=("1" "1")
	# # printf "\nPlotting median profile for ${nameout} scaling by sample\n"
	# # plotProfile -m manuscript/matrix/${nameout}.gz -out manuscript/plots/${nameout}_scaled_profile_median.pdf --samplesLabel ${labels[@]} "Genes" "TEs" --regionsLabel ${region_label[@]} --averageType median --yMin ${ymins[@]} --yMax ${ymaxs[@]} --colors ${colorsprofile[@]}
	
	# # rm -f manuscript/data/temp*
	# # rm -f manuscript/matrix/temp*
	# # rm -f manuscript/plots/${nameout}_scaled_profile.pdf
	# # rm -f manuscript/data/distal_${tissue}_*.txt
	
	# # printf "Computing mC matrix ${tissue}\n"
	# # computeMatrix scale-regions -R manuscript/data/${tissue}_sorted_scaled_type_enhancers.txt -S ${mcbw[@]} -bs 50 -b 5000 -a 5000 -m 1000 --sortRegions keep -p ${threads} -o manuscript/matrix/matrix_mC_${tissue}.gz
	# # printf "Plotting heatmap ${tissue}\n"
	# # plotHeatmap -m manuscript/matrix/matrix_mC_${tissue}.gz -out manuscript/plots/${nameout}_scaled_mC.pdf --sortRegions keep --samplesLabel ${mcsamplelab[@]} --regionsLabel ${region_label[@]} --colorMap 'Oranges' --missingDataColor 'grey' --interpolationMethod 'nearest' --yMin ${mcmins[@]} --yMax ${mcmaxs[@]} --zMin ${mcmins[@]} --zMax ${mcmaxs[@]} --xAxisLabel "enhancer"
	# # plotHeatmap -m manuscript/matrix/matrix_mC_${tissue}.gz -out manuscript/plots/${nameout}_scaled_heatmap_mC.pdf --sortRegions keep --samplesLabel ${mcsamplelab[@]} --regionsLabel ${region_label[@]} --colorMap 'Oranges' --missingDataColor 'grey' --interpolationMethod 'nearest' --yMin ${mcmins[@]} --yMax ${mcmaxs[@]} --zMin ${mcmins[@]} --zMax ${mcmaxs[@]} --xAxisLabel "enhancer" --whatToShow "heatmap only"
	# # plotProfile -m manuscript/matrix/matrix_mC_${tissue}.gz -out manuscript/plots/${nameout}_scaled_profile_mC_mean.pdf --samplesLabel ${mcsamplelab[@]} --regionsLabel ${region_label[@]} --numPlotsPerRow 3 --yMin ${mcmins[@]} --yMax ${mcmaxs[@]} --averageType mean --colors ${colorsprofile[@]}
	# # plotProfile -m manuscript/matrix/matrix_mC_${tissue}.gz -out manuscript/plots/${nameout}_scaled_profile_mC_median.pdf --samplesLabel ${mcsamplelab[@]} --regionsLabel ${region_label[@]} --numPlotsPerRow 3 --yMin ${mcmins[@]} --yMax ${mcmaxs[@]} --averageType median --colors ${colorsprofile[@]}
# done

###################################################################################################################
### To make browser shots at enhancer (Fig 4b)
###################################################################################################################

# printf "chr9\t128009612\t128041835\tinter2\n" >> manuscript/data/other_regions_B73.txt
# while read chr start stop peakID exp strand gid
# do	
	# type=$(grep -w "${peakID}" manuscript/data/B73_v5_cn_with*exp*enhancers.bed | awk -v FS="cn_|_enhancers" '{print $2}')
	# printf "${peakID} ${type}\n"
	# grep -w "${peakID}" ChIP/peaks/best_peaks_B73_cn_H3K27ac.bed | awk -v t=${type} '{print $1,$2-2000,$3+2000,t"_"$4}' >> manuscript/data/other_regions_B73.txt
# done < manuscript/data/top_loci_cn.bed

# for peakID in peak_56145b peak_61972 peak_47085d peak_19323c peak_19028b
# do
	# type=$(grep "${peakID}" manuscript/data/B73_v5_cn_with*exp*enhancers.bed | awk -v FS="cn_|_enhancers" '{print $2}')
	# printf "${peakID} ${type}\n"
	# grep "${peakID}" ChIP/peaks/best_peaks_B73_cn_H3K27ac.bed | awk -v t=${type} '{print $1,$2-2000,$3+2000,t"_"$4}' >> manuscript/data/other_regions_B73.txt
# done

# pathtobwrna="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/RNA/tracks"
# pathtobwshrna="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/shRNA/tracks"
# pathtobwchip="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/ChIP/tracks"

# for line in B73
# do
	# nameout="SupFig4_${line}"
	# printf "Preparing samplefile for ${nameout}\n"
	# case ${line} in
		# B73) 	analysisname="B73_on_B73_v5_all_genes"
				# ref="B73_v5";;
		# W22)	analysisname="W22_on_W22_v2_all_genes"
				# ref="W22_v2";;
		# NC350)	analysisname="NC350_on_NC350_NAM_all_genes"
				# ref="NC350_NAM";;
	# esac
		
	# ## To prepare the samplefile of tracks to plot browser with all tissues (input files)
	# rm -f manuscript/data/${nameout}_samplefile.txt
	# rm -f manuscript/data/temp_distal_${nameout}_samplefile.txt
	# for tissue in cn ears endosperm roots
	# do
		# sample="${line}_${tissue}"
		# case "${tissue}" in
			# cn)	backcolor="grey50";;
			# ears)	backcolor="blue";;
			# endosperm)	backcolor="red";;
			# roots)	backcolor="green";;
		# esac
		# for mark in RNAseq
		# do
			# trackcolor="#7716ab"
			# fillcolor1="#9067a6"
			# fillcolor2="#9067a6"
			# for strand in plus minus
			# do
				# case "${strand}" in
					# plus)	invert="no"
							# lims="0;1.1";;
					# minus)	invert="invert"
							# lims="-1.1;0";;
				# esac
				# group="${mark}_${strand}"
				# pathtobw="${pathtobwrna}/${sample}_${mark}_merged_${strand}.bw"
				# # printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\t${lims}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
				# printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
			# done
		# done
		# for mark in RAMPAGE
		# do
			# trackcolor="grey10"
			# fillcolor1="grey50"
			# fillcolor2="grey50"
			# for strand in plus minus
			# do
				# case "${strand}" in
					# plus)	invert="no"
							# lims="0;0.15";;
					# minus)	invert="invert"
							# lims="-0.15;0";;
				# esac
				# group="${mark}_${strand}"
				# pathtobw="${pathtobwrna}/${sample}_${mark}_merged_${strand}.bw"
				# printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
			# done
		# done
		# for mark in 24RNA
		# do
			# trackcolor="#C45100"
			# fillcolor1="#FB9E00"
			# fillcolor2="#FB9E00"
			# for strand in plus minus
			# do
				# case "${strand}" in
					# plus)	invert="no"
							# lims="0;1";;
					# minus)	invert="invert"
							# lims="-1;0";;
				# esac
				# group="${mark}_${strand}"
				# pathtobw="${pathtobwshrna}/${sample}_${mark}_merged_${strand}.bw"
				# printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
			# done
		# done
		# for mark in H3K27ac H3K4me1 H3K4me3
		# do
			# group="${mark}"
			# case "${mark}" in
				# H3K27ac)	trackcolor="#EE616E"
							# fillcolor1="#ebd8da"
							# fillcolor2="#ed939b"
							# lims="-0.15;0.65";;
				# H3K4me1)	trackcolor="#8D9BEE"
							# fillcolor1="#d5d8e8"
							# fillcolor2="#afb8ed"
							# lims="-0.18;0.25";;
				# H3K4me3)	trackcolor="#F1C062"
							# fillcolor1="#e8e1d3"
							# fillcolor2="#ebcd94"
							# lims="-0.2;0.8";;
			# esac
			# pathtobw="${pathtobwchip}/${sample}_${mark}_merged.bw"
			# printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\tno\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
		# done
	# done
	
	# ### To prepare the locus file
	# sort -k1,1 -k2,2n ChIP/peaks/best_peaks_${line}_ears_H3K27ac.bed > manuscript/data/temp_distal_peaks_${nameout}.bed
	# bedtools intersect -v -a manuscript/data/temp_distal_peaks_${nameout}.bed -b ChIP/tracks/${ref}_all_genes.bed > manuscript/data/temp2_distal_peaks_${nameout}.bed
	# sort -k1,1 -k2,2n ChIP/peaks/best_peaks_${line}_cn_H3K27ac.bed > manuscript/data/temp_distal_peaks_${nameout}_v2.bed
	# bedtools intersect -v -a manuscript/data/temp_distal_peaks_${nameout}_v2.bed -b ChIP/tracks/${ref}_all_genes.bed >> manuscript/data/temp2_distal_peaks_${nameout}.bed
	# printf "Preparing loci file for ${nameout}\n"
	# rm -f manuscript/data/${nameout}_samplefile.txt
	# awk '$1~"cn"' manuscript/data/temp_distal_${nameout}_samplefile.txt > manuscript/data/${nameout}_samplefile.txt
	# # cat manuscript/data/temp_distal_${nameout}_samplefile.txt > manuscript/data/${nameout}_samplefile.txt
	# rm -f manuscript/data/${nameout}_loci.txt
	# while read chr start stop ID
	# do
		# windowstart=$((start))
		# windowstop=$((stop))
		# printf "${chr}\t${windowstart}\t${windowstop}\n" > manuscript/data/temp_distal_window_${nameout}.bed
		# enstart=$(bedtools intersect -wa -a manuscript/data/temp2_distal_peaks_${nameout}.bed -b manuscript/data/temp_distal_window_${nameout}.bed | awk -v ORS="," '{print $2}')
		# enwidth=$(bedtools intersect -wa -a manuscript/data/temp2_distal_peaks_${nameout}.bed -b manuscript/data/temp_distal_window_${nameout}.bed | awk -v ORS="," '{print $3-$2}')
		# printf "${chr}:${windowstart}:${windowstop}\t${nameout}_${ID}\t${enstart}\t${enwidth}\n" >> manuscript/data/${nameout}_loci.txt
	# done < manuscript/data/other_regions_${line}.txt
	
	# printf "\nRunning plot_browser_input_files script on ${nameout} ${short}\n"
	# cd manuscript
	# qsub -sync y /grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/Plot_browser_input_files.sh -l ${line} -r data/${nameout}_loci.txt -f data/${nameout}_samplefile.txt
	# cd ..
		
# done

###########################################################################################################
### To check values of "super"-enhancers and STARR-seq data (in B73_v4) from Ricci et al. 2019 (Fig4h, Sup Fig 7)
###########################################################################################################

# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/B73_v4"
# # bedGraphToBigWig STARR_seq/GSE120304_STARR_B73_enhancer_activity_ratio.txt ${ref_dir}/chrom.sizes STARR_seq/STARR_seq_B73_leaf.bw

# cat manuscript/data/B73_v4_distal_*.bed | awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -i - > manuscript/data/all_B73_v4_distal_types_and_tissues_merged.txt
# awk -v OFS="\t" '$1~/^[0-9]/ {if ($2>2000) print $1,$2-2000,$3+2000; else print $1,0,$3+2000}' ChIP/tracks/B73_v4_all_genes.bed | sort -k1,1 -k2,2n | bedtools merge -i - > manuscript/data/all_B73_v4_genes_and_surroundings.txt
# cat manuscript/data/all_B73_v4_distal_types_and_tissues_merged.txt manuscript/data/all_B73_v4_genes_and_surroundings.txt | sort -k1,1 -k2,2n | bedtools merge -i - > manuscript/data/all_B73_v4_enhancers_and_genes_and_surroundings.txt
# bedtools complement -i manuscript/data/B73_v4_mappable_genome.bed -g ${ref_dir}/chrom.sizes > manuscript/data/B73_v4_unmappable_genome.bed

# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/B73_v4"

# printf "Tissue\tGroup\tID\tMax\tMedian\tMean\tSize\n" > manuscript/data/total_B73_v4_STARR_v3.txt
# for tissue in cn ears endosperm roots
# do
	# ## To add shuffled regions in mappable genome
	# # bedtools shuffle -noOverlapping -i manuscript/data/B73_v4_distal_${tissue}_without_k4_expdouble.bed -g ${ref_dir}/chrom.sizes -excl manuscript/data/B73_v4_unmappable_genome.bed | awk -v OFS="\t" -v t=${tissue} '{print $1,$2,$3,t"_shuffled_"NR,NR,"."}' > manuscript/data/B73_v4_distal_${tissue}_shuffle.bed
	
	# cat manuscript/data/B73_v4_distal_${tissue}_shuffle.bed > manuscript/data/B73_v4_distal_${tissue}_shuffle.bed
	# for group in with_k4_expdouble with_k4_expsingle with_k4_noexp without_k4_expdouble without_k4_expsingle without_k4_noexp shuffle
	# do
		# bedtools intersect -wao -a manuscript/data/B73_v4_distal_${tissue}_${group}.bed -b STARR_seq/GSE120304_STARR_B73_enhancer_activity_ratio.txt | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2,$3,$4,$10,$11}' | bedtools merge -c 4,5,5,5,6 -o distinct,max,median,sum,sum -i - | awk -v OFS="\t" -v t=${tissue} -v g=${group} '{print t,g,$4,$5,$6,$7/$8,$8}' >> manuscript/data/total_B73_v4_STARR_v3.txt
	# done
# done

# rm -f manuscript/data/total_B73_v4_STARR*.txt
# printf "Tissue\tGroup\tID\tMax\tMedian\tMean\tSize\n" > manuscript/data/total_B73_v4_STARR_v2.txt
# for tissue in cn ears endosperm roots
# do
	# for group in with_k4_expdouble with_k4_expsingle with_k4_noexp without_k4_expdouble without_k4_expsingle without_k4_noexp
	# do
		# printf "${tissue} enhancers ${group}\n"
		# bedtools intersect -wao -a manuscript/data/B73_v4_distal_${tissue}_${group}.bed -b STARR_seq/GSE120304_STARR_B73_enhancer_activity_ratio.txt | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2,$3,$4,$10,$11}' | bedtools merge -c 4,5,5,5,6 -o distinct,max,median,sum,sum -i - | awk -v OFS="\t" -v t=${tissue} -v g=${group} '{print t,g,$4,$5,$6,$7/$8,$8}' >> manuscript/data/total_B73_v4_STARR_v2.txt
	# done
# done

## Then plot with MaizeCode_R_STARR_seq.r

# Rscript --vanilla MaizeCode_R_STARR_seq.r manuscript/data/total_B73_v4_STARR_v2.txt STARRseq

########################################################################
### To do the alluvial plot of DEG conservation accross inbreds (Fig 5a)
########################################################################

# ### For B73_v5
# awk '$3=="zea_maysb73" {print $2,$4}' crossref/Ti11_orthologs.tsv > manuscript/data/TIL11_vs_B73_v5_crossref.txt

# ### For NC350
# awk '$3=="zea_maysnc350" {print $2,$4}' crossref/Ti11_orthologs.tsv > manuscript/data/TIL11_vs_NC350_NAM_crossref.txt

# ### For W22
# awk '$3=="zea_maysti11" {print $4,$2}' crossref/W22_orthologs.tsv > manuscript/data/TIL11_vs_W22_v2_crossref.txt

# rm -f manuscript/data/all_TIL11.txt
# rm -f manuscript/data/*other_DEGs.txt
# for tissue in ears cn pollen roots
# do
	# for deg in UP DOWN
	# do
		# printf "\nMerging ${tissue} ${deg} TIL11\n"
		# awk -v OFS="\t" -v t=${tissue} -v d=${deg} '{print $4,t,d,"TIL11",$4}' combined/DEG/only_TIL11_${tissue}_DEG_${deg}_Alluvial_on_TIL11_chrs_all_genes.bed | sort -u >> manuscript/data/all_TIL11.txt
	# done
	# for line in TIL11 B73 W22 NC350
	# do
		# case "${line}" in
			# TIL11)	ref="TIL11_chrs";;
			# B73)	ref="B73_v5";;
			# NC350)	ref="NC350_NAM";;
			# W22)	ref="W22_v2";;
		# esac
		# printf "\nMerging other DEGs in ${tissue} ${line}\n"
		# cat combined/DEG/DEG_Alluvial_on_${ref}_all_genes_${line}*${tissue}*.txt | awk -v OFS="\t" -v t=${tissue} '$1 != "Chr" {print $4,t}' | sort -u > manuscript/data/${line}_${tissue}_other_DEGs.txt
	# done
# done

# cat manuscript/data/TIL11_*_other_DEGs.txt > manuscript/data/TIL11_other_DEGs.txt
# awk -v OFS="\t" '{print $1,$2,$3}' manuscript/data/all_TIL11.txt > manuscript/data/unique_TIL11_IDs.txt
# cat manuscript/data/unique_TIL11_IDs.txt > manuscript/data/all_TIL11_IDs.txt

# while read GID tissue
# do
	# if grep -q "${GID}" manuscript/data/unique_TIL11_IDs.txt; then
		# continue
	# else
		# printf "${GID}\t${tissue}\tOther_DEG\n" >> manuscript/data/all_TIL11_IDs.txt
		# printf "${GID}\t${tissue}\tOther_DEG\tTIL11\t${GID}\n" >> manuscript/data/all_TIL11.txt
	# fi
# done < manuscript/data/TIL11_other_DEGs.txt

# printf "\nSearching for homologs\n"

# while read GID tissue deg
# do
	# for line in B73 W22 NC350
	# do
		# case "${line}" in
			# B73)	ref="B73_v5";;
			# NC350)	ref="NC350_NAM";;
			# W22)	ref="W22_v2";;
		# esac
		# if grep -q "${GID}" manuscript/data/TIL11_vs_${ref}_crossref.txt; then
			# grep "${GID}" manuscript/data/TIL11_vs_${ref}_crossref.txt | awk '{print $2}' > manuscript/data/tempb.txt
			# a="No"
			# while read tempid
			# do
				# if [[ ${deg} == "UP" ]] || [[ ${deg} == "DOWN" ]]; then
					# if [[ ${a} == "No" ]]; then
						# if grep -q "${tempid}" combined/DEG/only_${line}_${tissue}_DEG_${deg}_Alluvial_on_${ref}_all_genes.bed; then
							# expr="${deg}"
							# gid="${tempid}"
							# a="Yes"
						# elif grep -q "${tempid}" manuscript/data/${line}_${tissue}_other_DEGs.txt; then
							# expr="Other_DEG"
							# gid="${tempid}"
						# else
							# expr="homolog_noDEG_${deg}"
							# gid="${tempid}"
						# fi
					# fi
				# else 	
					# if [[ ${a} == "No" ]]; then
						# if grep -q "${tempid}" combined/DEG/only_${line}_${tissue}_DEG_UP_Alluvial_on_${ref}_all_genes.bed; then
							# expr="UP"
							# gid="${tempid}"
							# a="Yes"
						# elif grep -q "${tempid}" combined/DEG/only_${line}_${tissue}_DEG_DOWN_Alluvial_on_${ref}_all_genes.bed; then
							# expr="DOWN"
							# gid="${tempid}"
							# a="Yes"
						# elif grep -q "${tempid}" manuscript/data/${line}_${tissue}_other_DEGs.txt; then
							# expr="Other_DEG"
							# gid="${tempid}"
						# else
							# expr="homolog_noDEG_${deg}"
							# gid="${tempid}"
						# fi
					# fi
				# fi
			# done < manuscript/data/tempb.txt
		# else
			# gid="NA"
			# expr="no_homolog"
		# fi
		# printf "${GID}\t${tissue}\t${expr}\t${line}\t${gid}\n" >> manuscript/data/all_TIL11.txt
	# done
# done < manuscript/data/all_TIL11_IDs.txt

# printf "\nCreating alluvial plot for Fig5\n"
# Rscript --vanilla MaizeCode_R_Cross_alluvial.r manuscript/data/all_TIL11.txt


#########################################################################################################################
############### Conservation analysis (Fig 5b,c)
#########################################################################################################################

### Intersect enhancers with conserved regions from Armin and CNS called by Conservatory2

# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/B73_v5"

# rm -f manuscript/data/temp_table_enhancers_for_conservation.bed
# for tissue in cn ears endosperm roots
# do
	# for type in with_k4_expdouble with_k4_expsingle with_k4_noexp without_k4_expdouble without_k4_expsingle without_k4_noexp
	# do
		# awk -v OFS="\t" -v t=${tissue} -v g=${type} '{print $1,$2,$3,$4,t,g}' manuscript/data/B73_v5_${tissue}_${type}_enhancers.bed >> manuscript/data/temp_table_enhancers_for_conservation.bed
	# done
	# bedtools shuffle -noOverlapping -i manuscript/data/B73_v5_${tissue}_without_k4_expdouble_enhancers.bed -g  manuscript/data/B73_v5_chrom_sizes.txt -excl manuscript/data/all_B73_v5_distal_types_and_tissues_merged.txt | awk -v OFS="\t" -v t=${tissue} '{print $1,$2,$3,t"_shuffled_"NR,t,"shuffled"}' >> manuscript/data/temp_table_enhancers_for_conservation.bed
# done

# sort -k1,1 -k2,2n manuscript/data/temp_table_enhancers_for_conservation.bed > manuscript/data/all_enhancers_for_conservation.bed

# printf "ID\tTissue\tType\tNbConserved\n" > manuscript/data/table_B73_v5_enhancers_merged_with_conserv_raw.txt
# bedtools intersect -c -a manuscript/data/all_enhancers_for_conservation.bed -b Armin_conservation_pan_andropogoneae/B73_v5_conserved_regions_raw.bed | awk -v OFS="\t" '{print $4,$5,$6,$7}' >> manuscript/data/table_B73_v5_enhancers_merged_with_conserv_raw.txt

# bedtools merge -i Anat_ConservatoryCNS/conservatory_CNS.bed > Anat_ConservatoryCNS/conservatory_CNS_merged.bed
# printf "ID\tTissue\tType\tNbConserved\n" > manuscript/data/table_B73_v5_enhancers_with_CNS_merged.txt
# bedtools intersect -c -a manuscript/data/all_enhancers_for_conservation.bed -b Anat_ConservatoryCNS/conservatory_CNS_merged.bed | awk -v OFS="\t" '{print $4,$5,$6,$7}' >> manuscript/data/table_B73_v5_enhancers_with_CNS_merged.txt

# ### Intersect enhancers with CNS called by Conservatory2 by percentage coverage

# rm -f manuscript/data/table_B73_v5_enhancers_with_CNS_v2.txt
# printf "ID\tTissue\tType\tConserv_percentage\n" > manuscript/data/table_B73_v5_enhancers_with_CNS_v2.txt
# for tissue in cn ears endosperm roots
# do
	# for type in with_k4_expdouble with_k4_expsingle with_k4_noexp without_k4_expdouble without_k4_expsingle without_k4_noexp shuffled
	# do
		# if [[ "${type}" == shuffle ]]; then	
			# sort -k1,1 -k2,2n manuscript/data/B73_v5_${tissue}_${type}_enhancers.bed > manuscript/data/temp_B73_v5_${tissue}_${type}_enhancers.bed
		# else
			# awk -v OFS="\t" '{print $1,$2,$3,$4}' manuscript/data/B73_v5_${tissue}_${type}_enhancers.bed | sort -k1,1 -k2,2n > manuscript/data/temp_B73_v5_${tissue}_${type}_enhancers.bed
		# fi
		# bedtools intersect -wao -a manuscript/data/temp_B73_v5_${tissue}_${type}_enhancers.bed -b Anat_ConservatoryCNS/conservatory_CNS_merged.bed | awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' > manuscript/data/temp_table_${tissue}_${type}_with_CNS_v2.txt
		# bedtools merge -i manuscript/data/temp_table_${tissue}_${type}_with_CNS_v2.txt -c 4,5 -o distinct,sum | awk -v OFS="\t" '{print $1,$2,$3,$4,$5}' > manuscript/data/table_${tissue}_${type}_with_CNS_v2.txt
		# awk -v OFS="\t" -v t=${tissue} -v g=${type} '{print $4,t,g,$5/($3-$2+1)*100}' manuscript/data/table_${tissue}_${type}_with_CNS_v2.txt >> manuscript/data/table_B73_v5_enhancers_with_CNS_v2.txt
	# done	
# done

# rm -f manuscript/data/temp_*

# ## Then plot with 

# Rscript --vanilla - <<-'EOF'
	# #!/usr/bin/env Rscript
	# library(dplyr)
	# library(tidyr)
	# library(ggplot2)
	# library(ggpubr)

	# tab<-read.delim("manuscript/data/table_B73_v5_enhancers_merged_with_conserv_raw.txt", header = TRUE)
	# tab$Tissue<-as.factor(tab$Tissue)
	# tab$Type<-factor(tab$Type, levels=c("with_k4_expdouble", "with_k4_expsingle", "with_k4_noexp",
          # "without_k4_expdouble", "without_k4_expsingle", "without_k4_noexp","shuffled"))

	# tab2<-mutate(tab, Conservation=ifelse(NbConserved==0, "0",
                                      # ifelse(NbConserved<5, "1-4",
                                             # ifelse(NbConserved<10, "5-9", "10+"))))
	# tab2$Conservation<-factor(tab2$Conservation, levels=c("0","1-4","5-9", "10+"))
	# tab3<-tab2 %>% group_by(Tissue, Type) %>%
		# summarize(Tot=n())

	# tabtot<-merge(tab2,tab3, by=c("Tissue","Type"))
	
	# plot1<-ggplot(tab2, aes(Type, fill=Conservation)) +
		# geom_bar(position="fill") +
		# scale_fill_manual(values = c("0"="grey80","1-4"="gold","5-9"="orange","10+"="red")) +
		# facet_grid(~Tissue) +
		# geom_text(data=tabtot, aes(x=Type, label=Tot), y=0.9, size=2, color="black", angle=45) +
		# theme(axis.text.x = element_text(angle = 45, hjust = 1),
				# axis.text.y = element_text(size=6),
				# axis.title = element_blank(),
				# panel.grid=element_blank(),
				# panel.background=element_rect(fill=NA, color="black"),
				# strip.background=element_blank(),
				# legend.title=element_text(size=6),
				# legend.text=element_text(size=6),
				# legend.key.size=unit(3,"mm")) +
		# labs(fill="Number of conserved regions")

	# pdf("manuscript/plots/Fig5b_Conserved_regions_B73_enhancers_Armin.pdf", height = 5, width=8)
	# print(plot1)
	# dev.off()
	
	# # ###
	
	# tab<-read.delim("manuscript/data/table_B73_v5_enhancers_with_CNS_merged.txt", header = TRUE)
	# tab$Tissue<-as.factor(tab$Tissue)
	# tab$Type<-factor(tab$Type, levels=c("with_k4_expdouble", "with_k4_expsingle", "with_k4_noexp",
          # "without_k4_expdouble", "without_k4_expsingle", "without_k4_noexp","shuffled"))

	# tab2<-mutate(tab, Conservation=ifelse(NbConserved==0, "0",
                                      # ifelse(NbConserved<5, "1-4",
                                             # ifelse(NbConserved<10, "5-9", "10+"))))
	# tab2$Conservation<-factor(tab2$Conservation, levels=c("0","1-4","5-9", "10+"))
	# tab3<-tab2 %>% group_by(Tissue, Type) %>%
		# summarize(Tot=n())

	# tabtot<-merge(tab2,tab3, by=c("Tissue","Type"))
	
	# plot1<-ggplot(tab2, aes(Type, fill=Conservation)) +
		# geom_bar(position="fill") +
		# scale_fill_manual(values = c("0"="grey80","1-4"="gold","5-9"="orange","10+"="red")) +
		# facet_grid(~Tissue) +
		# geom_text(data=tabtot, aes(x=Type, label=Tot), y=0.9, size=2, color="black", angle=45) +
		# theme(axis.text.x = element_text(angle = 45, hjust = 1),
				# axis.text.y = element_text(size=6),
				# axis.title = element_blank(),
				# panel.grid=element_blank(),
				# panel.background=element_rect(fill=NA, color="black"),
				# strip.background=element_blank(),
				# legend.title=element_text(size=6),
				# legend.text=element_text(size=6),
				# legend.key.size=unit(3,"mm")) +
		# labs(fill="Number of conserved regions")

	# pdf("manuscript/plots/Fig5c_Conserved_regions_B73_enhancers_CNS_merged_Anat.pdf", height = 5, width=8)
	# print(plot1)
	# dev.off()
	
# EOF

############################################################################################################################
###### To make plot of sRNA sizes (Sup Fig 3)
############################################################################################################################

# rm -f manuscript/data/total_table_shRNA_sizes.txt
# printf "Line\tTissue\tRep\tSize\tNb\n" > manuscript/data/total_table_shRNA_sizes.txt
# for line in B73 W22 NC350 TIL11
# do
	# for tissue in cn ears endosperm pollen roots
	# do
		# for rep in Rep1 Rep2
		# do
			# if [ -s shRNA/reports/sizes_mapped_${line}_${tissue}_shRNA_${rep}.txt ]; then
				# awk -v OFS="\t" -v l=${line} -v t=${tissue} -v r=${rep} '{print l,t,r,$3,$4}' shRNA/reports/sizes_mapped_${line}_${tissue}_shRNA_${rep}.txt >> manuscript/data/total_table_shRNA_sizes.txt
			# fi
		# done
	# done
# done

# Rscript --vanilla - <<-'EOF'
	# #!/usr/bin/env Rscript

	# library(dplyr)
	# library(tidyr)
	# library(ggplot2)
	# library(cowplot)
	
	# ori<-read.delim("manuscript/data/total_table_shRNA_sizes.txt", header=TRUE) %>% 
		# mutate(Length=ifelse(Size<=17,"17-",
                       # ifelse(Size>=76,"76-150",
						# ifelse(Size>=31,"31-75",Size)))) %>%
		# group_by(Line, Tissue, Rep) %>%
		# summarize(Total=sum(Nb))

	# cpm<-read.delim("manuscript/data/total_table_shRNA_sizes.txt", header=TRUE) %>%
		# mutate(Length=ifelse(Size<=17,"17-",
                       # ifelse(Size>=76,"76-150",
					    # ifelse(Size>=31,"31-75",Size)))) %>%
		# group_by(Line, Tissue, Rep, Length) %>%
		# summarize(Inter=sum(Nb)) %>%
		# merge(ori, by=c("Line","Tissue","Rep")) %>%
		# rowwise() %>%
		# mutate(CPM=Inter/Total*1000000) %>%
		# ungroup() %>%
		# select(-Inter,-Total) %>%
		# mutate(Color=ifelse(Length==21,"red",
						# ifelse(Length==22,"blue",
                             # ifelse(Length==24,"orange","white")))) %>%
		# spread(key=Rep, value=CPM, fill = 0) %>%
		# mutate(Mean=((Rep1+Rep2)/2), SD=abs(Mean-Rep1))

	# partial<-filter(cpm, Length %in% seq(20,25))
	
	# partial$Line<-factor(partial$Line, c("B73","NC350","W22","TIL11"))
	# partial$Color<-as.factor(partial$Color)
	
	# plottissue<-function(tissue) {
		# tab<-filter(partial, Tissue==tissue)
  
		# plot<-ggplot(tab, aes(Length,Mean,fill=Color)) +
			# geom_bar(stat="identity", position="stack", color="black", linewidth=0.01, show.legend = FALSE) +
			# geom_errorbar(aes(Length, ymin=Mean-SD, ymax=Mean+SD), colour="black", linewidth=0.5, width=0.5) +
			# facet_wrap(~Line, nrow=1) +
			# scale_fill_manual(values = c("red"="red","blue"="blue","orange"="orange","white"="white")) +
			# scale_y_continuous(labels=function(x)x/10000, name = expression(CPM~(x10^4))) +
			# theme_bw() +
			# theme(axis.text.x = element_text(angle=45, vjust=0.5),
				# panel.grid.minor.y = element_blank(),
				# axis.title.x = element_blank()) +
			# labs(title=paste0(tissue))
  
		# print(plot)
	# }

	# p1<-plottissue("cn")
	# p2<-plottissue("ears")
	# p3<-plottissue("pollen")
	# p4<-plottissue("endosperm")
	# p5<-plottissue("roots")

	# p4b<-plot_grid(p4,NULL, rel_widths = c(3.15,1))
	# p5b<-plot_grid(p5,NULL, rel_widths = c(3.15,1))
	
	# plot5<-plot_grid(p1,p2,p3,p4b,p5b,nrow=5)
	
	# pdf(paste0("manuscript/plots/SupFig_shRNA_sizes_stats_short_scale.pdf"), height=10, width=12)
	# print(plot5)
	# dev.off()
	
# EOF


###########################################################################
### To make browser shots at b1 color gene (Sup Fig 5)
###########################################################################

# rm -f manuscript/data/color_genes_*.txt
# printf "b1\tcolored_plant1\tZm00001eb074320\tZm00004b006646\tZm00036ab076510\tcolored plant 1 (b1)\n" > manuscript/data/color_genes_all.txt

# # for line in B73 W22
# # do
	# pathtobwrna="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/RNA/tracks"
	# pathtobwshrna="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/shRNA/tracks"
	# pathtobwchip="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/ChIP/tracks"
	# pathtobwtf="/grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/TF/tracks"
	# nameout="color_genes_${line}"
	# ## To prepare the samplefile of tracks to plot browser with all tissues (input files)
	# rm -f manuscript/data/${nameout}_samplefile.txt
	# rm -f manuscript/data/temp_distal_${nameout}_samplefile.txt
	# rm -f manuscript/data/${nameout}_loci.txt
	# for tissue in cn ears
	# do
		# sample="${line}_${tissue}"
		# case "${tissue}" in
			# cn)	backcolor="grey50";;
			# ears)	backcolor="blue";;
			# endosperm)	backcolor="red";;
			# roots)	backcolor="green";;
		# esac
		# for mark in RNAseq
		# do
			# trackcolor="#7716ab"
			# fillcolor1="#9067a6"
			# fillcolor2="#9067a6"
			# for strand in plus minus
			# do
				# case "${strand}" in
					# plus)	invert="no"
							# lims="0;1.1";; ## for b1
					# minus)	invert="invert"
							# lims="-1.1;0";;	## for b1
				# esac
				# group="${mark}_${strand}"
				# pathtobw="${pathtobwrna}/${sample}_${mark}_merged_${strand}.bw"
				# # printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\t${lims}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
				# printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
			# done
		# done
		# for mark in RAMPAGE
		# do
			# trackcolor="grey10"
			# fillcolor1="grey50"
			# fillcolor2="grey50"
			# for strand in plus minus
			# do
				# case "${strand}" in
					# plus)	invert="no"
							# lims="0;0.15";;
					# minus)	invert="invert"
							# lims="-0.15;0";;
				# esac
				# group="${mark}_${strand}"
				# pathtobw="${pathtobwrna}/${sample}_${mark}_merged_${strand}.bw"
				# # printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\t${lims}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
				# printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
			# done
		# done
		# for mark in 24RNA shortRNA
		# do
			# case "${mark}" in
				# 24RNA)	trackcolor="orange"
						# fillcolor1="orange"
						# fillcolor2="orange";;
				# shortRNA)	trackcolor="brown"
						# fillcolor1="brown"
						# fillcolor2="brown";;			
			# esac
			# for strand in plus minus
			# do
				# case "${strand}" in
					# plus)	invert="no"
							# lims="0;0.15";;
					# minus)	invert="invert"
							# lims="-0.15;0";;
				# esac
				# group="${mark}_${strand}"
				# pathtobw="${pathtobwshrna}/${sample}_${mark}_merged_${strand}.bw"
				# # printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\t${lims}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
				# printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${invert}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
			# done
		# done
		# for mark in H3K27ac H3K4me1 H3K4me3
		# do
			# group="${mark}"
			# case "${mark}" in
				# H3K27ac)	trackcolor="#EE616E"
							# fillcolor1="#ebd8da"
							# fillcolor2="#ed939b"
							# lims="-0.15;0.65";;
				# H3K4me1)	trackcolor="#8D9BEE"
							# fillcolor1="#d5d8e8"
							# fillcolor2="#afb8ed"
							# lims="-0.18;0.25";;
				# H3K4me3)	trackcolor="#F1C062"
							# fillcolor1="#e8e1d3"
							# fillcolor2="#ebcd94"
							# lims="-0.2;0.8";;
			# esac
			# pathtobw="${pathtobwchip}/${sample}_${mark}_merged.bw"
			# # printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\tno\t${lims}\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
			# printf "${sample}\t${group}\t${pathtobw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\tno\n" >> manuscript/data/temp_distal_${nameout}_samplefile.txt
		# done
	# done
	# while read short long GID WID NID name
	# do
		# printf "Preparing samplefile for ${nameout}\n"
		# case ${line} in
			# B73) 	analysisname="B73_on_B73_v5_all_genes"
					# ref="B73_v5"
					# id="${GID}";;
			# W22)	analysisname="W22_on_W22_v2_all_genes"
					# ref="W22_v2"
					# id="${WID}";;
		# esac
		# if [[ ${id} == "NA" ]]; then
			# continue
		# else
			# ### To prepare the locus file
			# sort -k1,1 -k2,2n ChIP/peaks/best_peaks_${line}_ears_H3K27ac.bed > manuscript/data/temp_distal_peaks_${nameout}.bed
			# bedtools intersect -v -a manuscript/data/temp_distal_peaks_${nameout}.bed -b ChIP/tracks/${ref}_all_genes.bed > manuscript/data/temp2_distal_peaks_${nameout}.bed
			# sort -k1,1 -k2,2n ChIP/peaks/best_peaks_${line}_cn_H3K27ac.bed > manuscript/data/temp_distal_peaks_${nameout}_v2.bed
			# bedtools intersect -v -a manuscript/data/temp_distal_peaks_${nameout}_v2.bed -b ChIP/tracks/${ref}_all_genes.bed > manuscript/data/temp2_distal_peaks_${nameout}_v2.bed
			# printf "Preparing loci file for ${nameout} ${short}\n"
			# rm -f manuscript/data/${nameout}_samplefile.txt
			# chr=$(grep -w ${id} ChIP/tracks/${ref}_all_genes.bed | awk '{print $1}')
			# # awk '$1~"ears" || $1~"cn" || $1~"endosperm"' manuscript/data/temp_distal_${nameout}_samplefile.txt > manuscript/data/${nameout}_samplefile.txt
			# cat manuscript/data/temp_distal_${nameout}_samplefile.txt > manuscript/data/${nameout}_samplefile.txt			
			# min=40000
			# max=5000
			# bins=10
			# genestart=$(grep -w ${id} ChIP/tracks/${ref}_all_genes.bed | awk '{print $2}')
			# genestop=$(grep -w ${id} ChIP/tracks/${ref}_all_genes.bed | awk '{print $3}')
			# genewidth=$(grep -w ${id} ChIP/tracks/${ref}_all_genes.bed | awk '{print $3-$2}')
			# windowstart=$((genestart-min))
			# windowstop=$((genestop+max))
			# printf "${chr}\t${windowstart}\t${windowstop}\n" > manuscript/data/temp_distal_window_${nameout}.bed
			# enstart=$(bedtools intersect -wa -a manuscript/data/temp2_distal_peaks_${nameout}.bed -b manuscript/data/temp_distal_window_${nameout}.bed | awk -v ORS="," '{print $2}')
			# enwidth=$(bedtools intersect -wa -a manuscript/data/temp2_distal_peaks_${nameout}.bed -b manuscript/data/temp_distal_window_${nameout}.bed | awk -v ORS="," '{print $3-$2}')
			# enstart2=$(bedtools intersect -wa -a manuscript/data/temp2_distal_peaks_${nameout}_v2.bed -b manuscript/data/temp_distal_window_${nameout}.bed | awk -v ORS="," '{print $2}')
			# enwidth2=$(bedtools intersect -wa -a manuscript/data/temp2_distal_peaks_${nameout}_v2.bed -b manuscript/data/temp_distal_window_${nameout}.bed | awk -v ORS="," '{print $3-$2}')
			# insertstart=110019481
			# insertwidth=10
			# printf "${chr}:${windowstart}:${windowstop}\t${nameout}_${short}_${id}\t${bins}\t${genestart},${insertstart}\t${genewidth},${insertwidth}\n" >> manuscript/data/${nameout}_loci.txt
			# # if [[ "${enstart}" != "" ]]; then
				# # printf "${chr}:${windowstart}:${windowstop}\t${nameout}_${short}_${id}\t${bins}\t${genestart},${enstart}\t${genewidth},${enwidth}\n" >> manuscript/data/${nameout}_loci.txt
			# # else
				# # printf "${chr}:${windowstart}:${windowstop}\t${nameout}_${short}_${id}\t${bins}\t${genestart},${enstart2}\t${genewidth},${enwidth2}\n" >> manuscript/data/${nameout}_loci.txt
			# # fi
		# fi
	# # done < manuscript/data/color_genes_all.txt

	# printf "\nRunning plot_browser_input_files script on ${nameout}\n"
	# cd manuscript/
	# qsub -sync y /grid/martienssen/home/jcahn/nlsas/projects/MaizeCode/Plot_browser_input_files.sh -l ${line} -r data/${nameout}_loci.txt -f data/${nameout}_samplefile.txt
	# cd ..
# done



#####################################################################################################
### to look at overlap between distal enhancers (Sup Fig  6c)
#####################################################################################################

# for line in B73 W22 NC350
# do
	# case "${line}" in
		# B73)	ref="B73_v5";;
		# W22)	ref="W22_v2";;
		# NC350)	ref="NC350_NAM";;
	# esac
	# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays/${ref}"

	# rm -f manuscript/data/temp_distal_without_expdouble.txt
	# for tissue in cn ears endosperm roots
	# do
		# awk -v OFS="\t" -v t=${tissue} '{print $1,$2,$3,t,$4}' manuscript/data/${ref}_${tissue}_without_k4_expdouble_enhancers.bed >> manuscript/data/temp_distal_without_expdouble.txt
	# done
	# sort -k1,1 -k2,2n manuscript/data/temp_distal_without_expdouble.txt > manuscript/data/all_${ref}_without_expdouble.txt
	# rm -f manuscript/data/temp_distal_without_expdouble.txt

	# bedtools merge -i manuscript/data/all_${ref}_without_expdouble.txt -c 4 -o distinct | bedtools sort -g ${ref_dir}/chrom.sizes | awk -v OFS="\t" '{print $1,$2,$3,"enhancer_"NR,$4}' > manuscript/data/all_${ref}_distal_without_expdouble_merged.txt

	# rm -f manuscript/data/temp_col_ALL_merged_*.txt
	# for tissue in cn ears endosperm roots
	# do
		# printf "${tissue}\n" > manuscript/data/temp_col_ALL_merged_${tissue}.txt
		# awk -v OFS="\t" -v t=${tissue} '{if ($5 ~ t) print "1"; else print "0"}' manuscript/data/all_${ref}_distal_without_expdouble_merged.txt >> manuscript/data/temp_col_ALL_merged_${tissue}.txt
	# done

	# awk -v OFS="\t" 'BEGIN {print "ID"} {print $4}' manuscript/data/all_${ref}_distal_without_expdouble_merged.txt > manuscript/data/temp_col_ALL_merged_AAA.txt
	# paste manuscript/data/temp_col_ALL_merged_*.txt | uniq > manuscript/data/matrix_upset_without_expdouble.txt
	# rm -f manuscript/data/temp_col_ALL_merged_*.txt

	# #### To make an Upset plot
	# printf "\nCreating Upset plot for without_expdouble_${ref}\n"
	# Rscript --vanilla R_Upset.r manuscript/data/matrix_upset_without_expdouble.txt manuscript/plots/Fig6c_upset_without_expdouble_${ref}
	
# done

############################################################################################################################
### To look at overlap between distal enhancer types, OCRS and loops (Sup Fig 8a) and the gene expression of genes in loops
############################################################################################################################

### To look at OCRs in super-enhancers (SupFig 8a)

# for tissue in cn ears endosperm roots
# do	
	# for group in with_k4_expdouble with_k4_expsingle with_k4_noexp without_k4_expdouble without_k4_expsingle without_k4_noexp shuffled
	# do
		# awk -v OFS="\t" -v t=${tissue} -v g=${group} '{print $1,$2,$3,$4,t,g}' manuscript/data/B73_v4_distal_${tissue}_${group}.bed >> manuscript/data/temp_total_table_B73_v4_enhancers.bed
	# done
# done

# sort -k1,1 -k2,2n manuscript/data/temp_total_table_B73_v4_enhancers.bed > manuscript/data/total_table_B73_v4_enhancers.bed 
# rm -f manuscript/data/temp_total_table_B73_v4_enhancers.bed

# printf "ID\tTissue\tType\tOCRs\n" > manuscript/data/table_B73_v4_enhancers_with_Sun_OCRs_ears.txt
# bedtools intersect -c -a manuscript/data/total_table_B73_v4_enhancers.bed -b manuscript/data/Sun_OCRs_ears.bed | awk -v OFS="\t" '{print $4,$5,$6,$7}' >> manuscript/data/table_B73_v4_enhancers_with_Sun_OCRs_ears.txt

# printf "ID\tTissue\tType\tOCRs\n" > manuscript/data/table_B73_v4_enhancers_with_Sun_OCRs_tassel.txt
# bedtools intersect -c -a manuscript/data/total_table_B73_v4_enhancers.bed -b manuscript/data/Sun_OCRs_tassel.bed | awk -v OFS="\t" '{print $4,$5,$6,$7}' >> manuscript/data/table_B73_v4_enhancers_with_Sun_OCRs_tassel.txt

# printf "ID\tTissue\tType\tOCRs\n" > manuscript/data/table_B73_v4_enhancers_with_Schmitz_ATAC_ears.txt
# bedtools intersect -c -a manuscript/data/total_table_B73_v4_enhancers.bed -b manuscript/data/Schmitz_ATAC_ears.bed | awk -v OFS="\t" '{print $4,$5,$6,$7}' >> manuscript/data/table_B73_v4_enhancers_with_Schmitz_ATAC_ears.txt

# printf "ID\tTissue\tType\tOCRs\n" > manuscript/data/table_B73_v4_enhancers_with_Schmitz_ATAC_leaf.txt
# bedtools intersect -c -a manuscript/data/total_table_B73_v4_enhancers.bed -b manuscript/data/Schmitz_ATAC_leaf.bed | awk -v OFS="\t" '{print $4,$5,$6,$7}' >> manuscript/data/table_B73_v4_enhancers_with_Schmitz_ATAC_leaf.txt

# rm -f manuscript/data/table_B73_v4_enhancers_with_loops_GIDs.txt
# printf "ID\tTissue\tType\tGID\n" > manuscript/data/table_B73_v4_enhancers_with_loops_GIDs.txt
# for tissue in cn ears endosperm roots
# do
	# for type in with_k4_expdouble with_k4_expsingle with_k4_noexp without_k4_expdouble without_k4_expsingle without_k4_noexp shuffle
	# do
		# bedtools intersect -wao -a manuscript/data/B73_v4_distal_${tissue}_${type}.bed -b manuscript/data/Sun_loops_ears.bed | awk -v OFS="\t" -v t=${tissue} -v g=${type} '{print $4,t,g,$10}' | sort -u >> manuscript/data/table_B73_v4_enhancers_with_loops_GIDs.txt
	# done	
# done

# for tissue in cn ears endosperm roots
# do
	# nb=$(awk -v t=${tissue} '$2==t && $3=="without_k4_expdouble" && $4!="."' manuscript/data/table_B73_v4_enhancers_with_loops_GIDs.txt | wc -l)
	# shuf -n ${nb} manuscript/data/Sun_loops_ears.bed | awk -v OFS="\t" -v t=${tissue} '{print "Random_loops"NR,t,"random_loops",$4}' >> manuscript/data/table_B73_v4_enhancers_with_loops_GIDs.txt
# done

# rm -f manuscript/data/temp*

### Then plot with R

# Rscript --vanilla - <<-'EOF'
	# #!/usr/bin/env Rscript
	# library(dplyr)
	# library(tidyr)
	# library(ggplot2)
	# library(stringr)
	# library(ggrepel)
	# library(ggpubr)
	# library(data.table)
	# library(ggalluvial)

	# schmitz_leaf<-read.delim("Input/table_B73_v4_enhancers_with_Schmitz_ATAC_leaf.txt", header=TRUE) %>%
	#   rename(ricci_leaf=OCRs)
	# schmitz_ears<-read.delim("Input/table_B73_v4_enhancers_with_Schmitz_ATAC_ears.txt", header=TRUE) %>%
	#   rename(ricci_ears=OCRs)
	# sun_ears<-read.delim("Input/table_B73_v4_enhancers_with_Sun_OCRs_ears.txt", header=TRUE) %>%
	#   rename(sun_ears=OCRs)
	# sun_tassel<-read.delim("Input/table_B73_v4_enhancers_with_Sun_OCRs_tassel.txt", header=TRUE) %>%
	#   rename(sun_tassel=OCRs)

	# sun_loops<-read.delim("Input/table_B73_v4_enhancers_with_loops_GIDs.txt", header=TRUE) %>%
	#   mutate(loops=ifelse(GID==".",0,1)) %>%
	#   group_by(ID,Tissue,Type) %>%
	#   summarize(sun_loops=sum(loops))

	# all<-merge(schmitz_ears,schmitz_leaf, by=c("ID","Tissue","Type")) %>%
	#   merge(sun_ears, by=c("ID","Tissue","Type")) %>%
	#   merge(sun_tassel, by=c("ID","Tissue","Type")) %>%
	#   merge(sun_loops, by=c("ID","Tissue","Type")) %>%
	#   mutate(regions_ricci_ears=ifelse(ricci_ears>=5, "5+", ricci_ears),
	#          regions_ricci_leaf=ifelse(ricci_leaf>=5, "5+", ricci_leaf),
	#          regions_sun_ears=ifelse(sun_ears>=5, "5+", sun_ears),
	#          regions_sun_loops=ifelse(sun_loops>=5, "5+", sun_loops),
	#          regions_sun_tassel=ifelse(sun_tassel>=5, "5+", sun_tassel)) %>%
	#   select(-ricci_ears, -ricci_leaf, -sun_ears, -sun_loops, -sun_tassel) %>%
	#   rename(ricci_ears=regions_ricci_ears, ricci_leaf=regions_ricci_leaf, sun_ears=regions_sun_ears, sun_loops=regions_sun_loops, sun_tassel=regions_sun_tassel)

	# all$Tissue<-as.factor(all$Tissue)
	# all$Type<-factor(all$Type, levels=c("with_k4_expdouble", "with_k4_expsingle", "with_k4_noexp",
	#                                             "without_k4_expdouble", "without_k4_expsingle", "without_k4_noexp","shuffle"))

	# temp<-group_by(all, ricci_ears, ricci_leaf, sun_ears, sun_loops, sun_tassel, Tissue, Type) %>%
	#   summarize(Counts=n()) %>%
	#   tibble::rownames_to_column("Groups") %>%
	#   gather(key="OCR",value="Nb",ricci_ears, ricci_leaf, sun_ears, sun_tassel, sun_loops)

	# temp$OCR<-factor(temp$OCR, levels=c("sun_tassel", "sun_ears", "sun_loops", "ricci_ears","ricci_leaf"))

	# temp$Type<-factor(temp$Type, levels=c("with_k4_expdouble", "with_k4_expsingle", "with_k4_noexp",
        #                   "without_k4_expdouble", "without_k4_expsingle", "without_k4_noexp","shuffle"))
	# temp$Groups<-as.factor(temp$Groups)

	# plot.OCRs.v2<-function(type) {
	  # tab<-filter(temp, Type==type)
  
	  # plot<-ggplot(tab, aes(x=OCR, stratum=Nb, alluvium=Groups, fill=Nb, y=Counts)) +
	    # geom_flow(aes.flow = "backward") +
	    # geom_stratum(width=0.5, size=0.1) +
	    # facet_wrap(~Tissue, scales="free_y", ncol=1, strip.position = "right") +
	    # scale_fill_manual(values = c("0"="grey80","1"="gold","2"="orange","3"="red", "4"="maroon", "5+"="black")) +
	    # labs(title=paste0(type)) +
	    # scale_x_discrete(expand = c(0,0)) +
	    # theme(panel.grid.minor = element_blank(),
	          # panel.grid.major.x = element_blank(),
	          # panel.grid.major.y = element_blank(),
	          # panel.background = element_blank(),
        	  # strip.background = element_rect(fill="white", color="black"),
        	  # axis.ticks.x=element_blank(), 
        	  # axis.title=element_blank(),
        	  # plot.title = element_text(hjust = 0.5, size=15),
        	  # axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=1),
       	 	  # axis.text.y=element_text(size=10),
         	 # legend.position="none")
  
	  # print(plot)
	# }

	#### Sup fig 8a
	# pdf("Plots/Alluvial_enhancers_B73_v4_in_OCRs_v2.pdf", height=12, width=15)
	# ggarrange(plot.OCRs.v2("with_k4_expdouble"),plot.OCRs.v2("with_k4_expsingle"),plot.OCRs.v2("with_k4_noexp"),
        	  # plot.OCRs.v2("without_k4_expdouble"),plot.OCRs.v2("without_k4_expsingle"),plot.OCRs.v2("without_k4_noexp"), 
        	  # plot.OCRs.v2("shuffle"), nrow = 1)
	# dev.off()

 	####

  	# tableRPKM<-read.delim("Input/genes_B73_v4_rpkm_DEGs.txt", header=TRUE)

	# mergedt<-filter(sun_loops, Type!="shuffle", GID != ".") %>%
		# merge(tableRPKM, by=c("GID","Tissue"), all.x=TRUE)

	## To add random genes

	# mergedtot<-group_by(mergedt, Tissue, Type) %>%
		  # summarize(Total=n())
	# shuf<-filter(mergedtot, Type=="without_k4_expdouble")
	# for ( tissue in c("ears","cn","endosperm","roots") ) {
  		# temptab<-filter(tableRPKM, Tissue==tissue)
  		# nb<-filter(shuf, Tissue==tissue)$Total  
  		# shuftemp<-as.data.table(sapply(temptab, sample, nb)) %>%
    		# mutate(ID=paste0("Random_gene",row_number()), Type="Random_genes")
  		# mergedt<-rbind(mergedt, shuftemp)
	# }

	# mergedt$RPKM<-as.numeric(mergedt$RPKM)
	# mergedt$Tissue<-as.factor(mergedt$Tissue)
	# mergedt$Type<-factor(mergedt$Type, levels=c("with_k4_expdouble", "with_k4_expsingle", "with_k4_noexp",
                                            "without_k4_expdouble", "without_k4_expsingle", "without_k4_noexp","random_loops","Random_genes"))

	# mergedtot2<-group_by(mergedt, Tissue, Type) %>%
  		# summarize(Total=n())

	# comps<-list( c("without_k4_expdouble", "without_k4_expsingle"),
             	# c("without_k4_expdouble", "without_k4_noexp"),
              	# c("without_k4_expdouble", "random_loops"),
             	# c("without_k4_expdouble", "Random_genes"))

	# plot2<-ggplot(mergedt, aes(Type, log2(RPKM+0.1), fill=Type)) +
  		# geom_boxplot(position="dodge", notch = FALSE) +
  		# facet_grid(~Tissue) +
  		# geom_text(data=mergedtot2, aes(x=Type, label=Total), y=11, size=3, color="black", angle=45) +
  		# scale_fill_manual(values=c("with_k4_expdouble" = "#5546EE", "with_k4_expsingle" = "#33A9ED", "with_k4_noexp" = "#33CBED", 
                              	# "without_k4_expdouble" = "#ED4433", "without_k4_expsingle" = "#ED8A33", "without_k4_noexp" = "#EDB533",
                              	# "random_Genes" = "#A9A9A9","random_loops" = "grey90")) +
  		# theme_bw() +
  		# stat_compare_means(comparisons = comps, method = "t.test", label="p.format",
                     	# tip.length = c(0)) +
  		# theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        		# legend.position = "none")

	#### Sup fig 8c
	# pdf("Plots/Enhancers_B73_v4_loops_to_GID_expression.pdf", height=6, width=8)
	# print(plot2)
	# dev.off()
 
# EOF

#####################################################################
### To get percentage of K27ac peaks in loops vs OCRs from Sun 2020 (Table in Sup Fig 8b)

# awk -v OFS="\t" 'NR>1 {print $1,$2,$3,$11}' manuscript/data/Sun_2020_ears_loops.txt > manuscript/data/Sun_ears_left_anchors.bed
# awk -v OFS="\t" 'NR>1 {print $4,$5,$6,$11}' manuscript/data/Sun_2020_ears_loops.txt > manuscript/data/Sun_ears_right_anchors.bed

# cat manuscript/data/Sun_ears_left_anchors.bed manuscript/data/Sun_ears_right_anchors.bed | sort -k1,1 -k2,2n > manuscript/data/Sun_ears_all_anchors.bed
# for loop in Gene-gene Intergenic-gene Intergenic-intergenic
# do
	# grep "${loop}" manuscript/data/Sun_ears_all_anchors.bed > manuscript/data/Sun_ears_all_anchors_${loop}.bed
# done

# rm -f manuscript/data/table_enhancers_in_loop_types.txt

# printf "Enhancer_Type\tLoop_Type\tNb\n" > manuscript/data/table_enhancers_in_loop_types.txt

# for type in with_k4_expdouble with_k4_expsingle with_k4_noexp without_k4_expdouble without_k4_expsingle without_k4_noexp shuffled
# do
	# wc -l manuscript/data/B73_v4_distal_ears_${type}.bed | awk -v OFS="\t" -v t=${type} '{print t,"Total",$1}' >> manuscript/data/table_enhancers_in_loop_types.txt
	# for loop in Gene-gene Intergenic-gene Intergenic-intergenic
	# do
		# sort -k1,1 -k2,2n manuscript/data/B73_v4_distal_ears_${type}.bed | bedtools intersect -c -a - -b manuscript/data/Sun_ears_all_anchors_${loop}.bed | awk -v OFS="\t" -v g=${type} -v l=${loop} '{if ($7!=0) n+=1} END {print g,l,n}' >> manuscript/data/table_enhancers_in_loop_types.txt
	# done
# done

# for type in dOCR LoOCR
# do
	# grep "${type}" manuscript/data/Sun_OCRs_ears.bed | awk -v OFS="\t" -v t=${type} '{print $1,$2,$3,"ear_"t"_"NR,NR,"."}' > manuscript/data/Sun_OCRs_ears_prep_${type}.bed
	# wc -l manuscript/data/Sun_OCRs_ears_prep_${type}.bed | awk -v OFS="\t" -v t=${type} '{print t,"Total",$1}' >> manuscript/data/table_enhancers_in_loop_types.txt
	# for loop in Gene-gene Intergenic-gene Intergenic-intergenic
	# do
		# sort -k1,1 -k2,2n manuscript/data/Sun_OCRs_ears_prep_${type}.bed | bedtools intersect -c -a - -b manuscript/data/Sun_ears_all_anchors_${loop}.bed | awk -v OFS="\t" -v g=${type} -v l=${loop} '{if ($7!=0) n+=1} END {print g,l,n}' >> manuscript/data/table_enhancers_in_loop_types.txt
	# done
# done

# rm -f manuscript/data/B73_v4_local_ears.bed
# awk -v OFS="\t" '{print $1,$2,$3,$4,$5,"."}' combined/peaks/enhancers_genic_B73_v4_on_B73_v4_all_genes_ears.txt >> manuscript/data/B73_v4_local_ears.bed
# awk -v OFS="\t" '{print $1,$2,$3,$4,$5,"."}' combined/peaks/enhancers_promoter_B73_v4_on_B73_v4_all_genes_ears.txt >> manuscript/data/B73_v4_local_ears.bed
# awk -v OFS="\t" '{print $1,$2,$3,$4,$5,"."}' combined/peaks/enhancers_terminator_B73_v4_on_B73_v4_all_genes_ears.txt >> manuscript/data/B73_v4_local_ears.bed
# wc -l manuscript/data/B73_v4_local_ears.bed | awk -v OFS="\t" -v t="Local_K27ac_peaks" '{print t,"Total",$1}' >> manuscript/data/table_enhancers_in_loop_types.txt
# for loop in Gene-gene Intergenic-gene Intergenic-intergenic
# do
	# sort -k1,1 -k2,2n manuscript/data/B73_v4_local_ears.bed | bedtools intersect -c -a - -b manuscript/data/Sun_ears_all_anchors_${loop}.bed | awk -v OFS="\t" -v g="Local_K27ac_peaks" -v l=${loop} '{if ($7!=0) n+=1} END {print g,l,n}' >> manuscript/data/table_enhancers_in_loop_types.txt
# done


############################################################################################################################
### To generate table of all DEGs (Sup Data 1)
############################################################################################################################

# for line in B73 W22 NC350 TIL11
# do
	# print "Generating DEG table for ${line}\n"
	# printf "Chr\tStart\tStop\tGID\tlogFC\tStrand\tlogCPM\tPValue\tFDR\tSamples\tDEG\n" > manuscript/data/SupTable2_${line}.tab
	# for tissue1 in cn ears endosperm pollen roots
	# do
		# for tissue2 in cn ears endosperm pollen roots
		# do
			# if [ -e combined/DEG/DEG_${line}*${tissue1}*vs*${tissue2}* ]; then
				# awk -v OFS="\t" 'NR>1' combined/DEG/DEG_${line}*${tissue1}*vs*${tissue2}* >> manuscript/data/SupTable2_${line}.tab
			# fi
		# done
	# done
# done

############################################################################################################################
#### To generate mappable region of genomes
############################################################################################################################

# ref="B73_v5"
# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Zea_mays"
# fasta=${ref_dir}/${ref}/temp_shRNA_B73_v5.fa

# while read chr stop
# do
        # # if [[ ${chr} =~ ^[0-9] ]]; then
		# if [[ ${chr} =~ ^chr ]]; then
				# printf "\nGetting bed file of reads for chromosome ${chr}\n\n"
                # awk -v OFS="\t" -v c=${chr} -v s=${stop} 'BEGIN {for (i=1; i<=s-150; i=i+10) print c,i,i+150}' > manuscript/data/${ref}_all_reads_${chr}.bed
				# printf "\nGetting fasta file of reads for chromosome ${chr}\n\n"
				# seqtk subseq ${fasta} manuscript/data/${ref}_all_reads_${chr}.bed | gzip >> manuscript/data/${ref}_all_reads_${chr}.fa.gz
        # fi
# done < ${ref_dir}/${ref}/chrom.sizes

# zcat manuscript/data/${ref}_all_reads_*.fa.gz | gzip > manuscript/data/${ref}_all_reads.fa.gz

# printf "\nMapping with bowtie 2\n\n"
# bowtie2 -p ${threads} -f --end-to-end --met-file manuscript/matrix/bt2_${ref}_all_reads.txt -x ${ref_dir}/${ref}/${ref} -U manuscript/data/${ref}_all_reads.fa.gz -S manuscript/data/${ref}_all_reads.sam |& tee manuscript/matrix/mapping_${ref}_all_reads.txt
# printf "\nSamtools processing of sam/bam\n\n"
# samtools view -@ ${threads} -b -h -q 10 -F 256 -o manuscript/data/temp1_${ref}_all_reads.bam manuscript/data/${ref}_all_reads.sam
# samtools fixmate -@ ${threads} -m manuscript/data/temp1_${ref}_all_reads.bam manuscript/data/temp2_${ref}_all_reads.bam
# samtools sort -@ ${threads} -o manuscript/data/temp3_${ref}_all_reads.bam manuscript/data/temp2_${ref}_all_reads.bam
# samtools markdup -r -s -f manuscript/matrix/markdup_${ref}_all_reads.txt -@ ${threads} manuscript/data/temp3_${ref}_all_reads.bam manuscript/data/${ref}_all_reads.bam
# samtools index -@ ${threads} manuscript/data/${ref}_all_reads.bam
# printf "\nGetting coverage of bam\n\n"
# bamCoverage --binSize 1 --outFileFormat bedgraph --normalizeUsing None -b manuscript/data/${ref}_all_reads.bam -o manuscript/data/${ref}_all_reads.bg
# printf "\nGenerating final bed file\n\n"
# awk -v OFS="\t" '$4>0' manuscript/data/${ref}_all_reads.bg | bedtools merge -i - > manuscript/data/${ref}_mappable_genome.bed

# bedtools complement -i manuscript/data/${ref}_mappable_genome.bed -g ${ref_dir}/chrom.sizes > manuscript/data/${ref}_unmappable_genome.bed

# rm -f manuscript/data/${ref}_all_reads.sam
# rm -f manuscript/data/temp*


############################################################################################################################

