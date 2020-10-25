# maize-code

MaizeCode Pipeline Help


Step-by-Step pipeline

1) Make a MaizeCode folder somewhere
2) Check for packages requirements (see package_versions file for versions that work for sure) and install the ones you don’t have:
pigz; samtools; bowtie2; STAR; fastqc; cutadapt; bedtools; deeptools; macs2; idr; R (and R libraries: readr, ggplot2, UpSetR)
3) Organize your reference genome folders so that they each contain ONE fasta file and ONE gff3 file (has to have 'gene' in column 3 and exons must be linked by 'Parent' in column 9)
4) Copy all these scripts in ~/data/Scripts/. If you would rather have them in another folder, check in the MaizeCode.sh script and replace the location.
5) Make the samplefiles you want. An example of a samplefile is in the data folder (B73_endosperm_samplefile.txt) and a quick way to make them is at the bottom of the MaizeCode.sh file
6) Submit the MaizeCode.sh script, giving as argument -f the samplefile.txt of your choice and -p the path to your directory that contains the different genome directories.
7) Make the analysis_samplefiles you want. An example of an analysis samplefile is in the data folder (B73_endosperm_analysis_samplefile.txt) and a quick way to make them is at the bottom of the MaizeCode_ChIP_analysis.sh file
8) Submit the MaizeCode_analysis.sh script, giving as argument -f the analysis_samplefile.txt and -r the regions (bed file) to be plotted on. -s can be set if the analysis should stop after calling peaks and making bigwig files. It also stops there if the regionfile is missing. !!  DOES NOT WORK RIGHT NOW, NEED TO PROVIDE REGION FILE !!


Comments

- Check out the usage of each script before as it can give some information about little things that might need to be edited. Submitting each script without arguments (or followed by help) will return its usage.
- The RNA pipeline is not ready yet
- It should work for both Single-end or Paired-end data but the SE part has not been tested (and I might not have edited it well as I was changing the PE part). A non-issue for now since all the ChIP data is PE, but to keep in mind for potential future use.
- A lot of the steps are only happening in some cases that should be more or less self-explanatory.
- The whole pipeline creates a lot of report files and probably files that are not necessary to keep but for now I keep them like this for potential troubleshooting.
- For now I’ve used the MaizeCode.sh script for 8 samples at a time (all ChIPs from a specific Line x Tissue). It runs in a couple hours (depending on the size of the files). It processes each sample in parallel and asks for 20 threads for each so better to limit the samplefiles to this amount.
- Always process the Input samples with their corresponding ChIP in the MaizeCode.sh script. (It can also be done separately, before or after) but they need to be done for the MaizeCode_ChIP_analysis.sh script to run successfully.
- The MaizeCode_analysis part will have to be tweaked for each analysis, but a default use of it should give a first look at the data.
- These are preliminary version of the scripts!


Scripts description

- MaizeCode.sh
Create the different folders
Create indexes (either bt2 or star) and reference files (if needed)
Copy fastq files from their folder in seq/ to the fastq/ folder (if not already done)
Run an instance of MaizeCode_ChIP_sample.sh or MaizeCode_RNA_sample.sh for each sample

- MaizeCode_ChIP_analysis.sh
Run fastQC on the raw data
Trim adapters, low quality and small reads (<20bp) with cutadapt
Run fastQC on trimmed data
Map with bowtie2
Removes PCR duplicates with samtools
Gets some mapping stats (in reports/summary_mapping_stats.txt)

- MaizeCode_RNA_sample.sh
NOT DONE YET, but expectations are:
Run fastQC on the raw data
Trim adapters and low quality with cutadapt
Run fastQC on trimmed data
Map with STAR with different settings depending on the type of data (RNAseq, shRNA, RAMPAGE)

- MaizeCode_analysis.sh
Split the analysis_samplefile into ChIP or RNA samples
Processed each datatype individually and wait for them to finish (WAITING IS NOT WORKING AT THE MOMENT)
Combined the individual datatype analysis to make more plots (NOTHING THERE YET)

- MaizeCode_ChIP_analysis.sh
Take either each replicate independently or merge the bam files of the two replicates
Call peaks with macs2 (calls broad peaks for H3K4me1, and narrow peaks for H3K4me3 and H3K27ac)
Make bigwig files with deeptools (log2 FC vs Input, normalizing each file by CPM) 
Do IDR analysis of the replicates with idr (if possible)
Get intersection of peaks to create different groups (presence/absence of each mark) with bedtools merge/intersect
Get distance of peaks in each group to closest gene with bedtools closest
Create Upset plot to show overlap among the different samples (using MaizeCode_R_Upset.r script)

- MaizeCode_RNA_analysis.sh
NOT DONE YET, but expectations are:
Call differentially expressed genes between samples

- MaizeCode_R_Upset.r
Create an Upset plot of overlapping peaks and their presence in gene bodies (e.g. data/Upset_B73_endosperm.pdf)
