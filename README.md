# maize-code

MaizeCode Pipeline Help


Step-by-Step pipeline

1) Make a MaizeCode folder somewhere
2) Copy the MaizeCode scripts folder somewhere
3) Check for packages requirements (see package_versions file for versions that work for sure) and install the ones you don’t have:
pigz; samtools; bowtie2; STAR; fastqc; cutadapt; bedtools; deeptools; macs2; idr; R (and R libraries: readr, ggplot2, UpSetR)
4) Organize your reference genome folders so that they are all in the same main folder and that each contain ONE fasta file and ONE gff3 file (has to have 'gene' in column 3 and exons must be linked by 'Parent' in column 9)
5) Make the samplefiles you want. An example of a samplefile is in the data folder (B73_endosperm_samplefile.txt) and a quick way to make them is at the bottom of the MaizeCode.sh file
6) Submit the MaizeCode.sh script, giving as argument -f the samplefile.txt of your choice and -p the path to your directory that contains the different genome directories.
7) By default, it will proceed with the analysis. It can be set so that it does not proceed with the analysis at all (-s argument to the MaizeCode.sh submission) and it will not do the combined analysis if several different references are being used for mapping
8) If the analysis has not proceeded, or different samples want to be analyzed together (or merged samples are to be used for example), make the corresponding analysis_samplefile you want. An example of an analysis samplefile is in the data folder (B73_endosperm_analysis_samplefile.txt) and a quick way to make them is at the bottom of the MaizeCode_analysis.sh file
9) Submit the MaizeCode_analysis.sh script, giving as argument -f the analysis_samplefile.txt and -r the regions (bed file) to be plotted on. -s can be set if the analysis should stop after calling peaks and making bigwig files. It also stops there if the regionfile is missing.


Comments

- There is one wrapper script (MaizeCode.sh) that launches sub-scripts depending on what needs to be done.
- The MaizeCode_analysis script is called by default by the MaizeCode.sh script but can be used seperately for any additional analysis.
- Check out the usage of each script before as it can give some information about the details and requirement for the scripts. Submitting each script without arguments (or followed by -h) will return its usage.
- The RNA pipeline is not ready yet
- It should work for both Single-end or Paired-end data but the SE part has not been tested (and I might not have edited it well as I was changing the PE part). A non-issue for now since all the ChIP data is PE, but to keep in mind for potential future use.
- The whole pipeline creates a lot of report files and probably files that are not necessary to keep but for now I keep them like this for potential troubleshooting.
- For now I’ve used the MaizeCode.sh script for 8 samples at a time (all ChIPs from a specific Line x Tissue). It runs in a couple hours (depending on the size of the files). It processes each sample in parallel and asks for 20 threads for each so better to limit the samplefiles to this amount.
- Always process the Input samples with their corresponding ChIP in the MaizeCode.sh script. (It can also be done separately, before or after) but they need to be done for the MaizeCode_ChIP_analysis.sh script to run successfully.
- The MaizeCode_analysis part will have to be tweaked for each analysis, but a default use of it should give a first look at the data.
- These are preliminary version of the scripts!


Scripts description

- MaizeCode.sh - wrapper script for the whole pipeline
Creates the different folders
Starts the MaizeCode_check_environment.sh script for each environment (datatype * reference) that needs to be created
Copies fastq files from their folder in seq/ to the fastq/ folder (if not already done)
Runs an instance of MaizeCode_ChIP_sample.sh or MaizeCode_RNA_sample.sh for each sample
Waits for the samples to be mapped
Launches the MaizeCode_analysis.sh script if the -s argument (that stops after mapping) has not been given
If mapping all the files to the same reference genome, it will use the all_genes.bed file created by the check_environment script
Otherwise, it will not performed the combined analysis (it would not provide a regionfile)

- MaizeCode_check_environment.sh
Checks if here is ONE fasta and ONE gff3 file in the reference folder (and unzip them if required)
Makes a chrom.sizes file if not there (can be useful down the line for bedGraphtoBigWig for example)
Makes a all_genes.bed file if not there (will be used for analysis/plots)
Create the template for the stat files ([ChIP|RNA]/reports/summary_mapping_stats.txt)
Makes the bowtie2 or STAR indexes (for ChIP and RNA, respectively) if not already there

- MaizeCode_ChIP_sample.sh
Runs fastQC on the raw data
Trims adapters, low quality and small reads (<20bp) with cutadapt
Runs fastQC on trimmed data
Maps with bowtie2
Removes PCR duplicates with samtools
Gets some mapping stats (in ChIP/reports/summary_mapping_stats.txt)

- MaizeCode_RNA_sample.sh - NOT DONE YET, but expectations are:
Runs fastQC on the raw data
Trims adapters and low quality with cutadapt
Runs fastQC on trimmed data
Maps with STAR with different settings depending on the type of data (RNAseq, shRNA, RAMPAGE)
Gets some mapping stats (in RNA/reports/summary_mapping_stats.txt)

- MaizeCode_analysis.sh - wrapper script for the analysis pipeline
If new files are to be analyzed, send each group of sample of the same datatype to MaizeCode_ChIP_analysis.sh or MaizeCode_RNA_analysis.sh
Launches the MaizeCode_combined_analysis.sh if the -s option is not set (or no regionfile is given), 

- MaizeCode_ChIP_analysis.sh
Takes either each replicate independently or merge the bam files of the two replicates
Calls peaks with macs2 (calls broad peaks for H3K4me1, and narrow peaks for H3K4me3 and H3K27ac)
Makes bigwig files with deeptools (log2 FC vs Input, normalizing each file by CPM) 

- MaizeCode_RNA_analysis.sh - NOT DONE YET, but expectations are:
Makes bigwig files with deeptools (normalized by CPM)

- MaizeCode_combined_analysis.sh
Splits the samplefile into ChIP and RNA samples
For ChIP samples:
Does IDR analysis of both replicates with idr (if both are present in the samplefile)
Makes a single file, merging all peaks from all samples with bedtools merge
Gets distance of each peak to the closest region from the regionfile with bedtools closest (default all genes annotated in the reference)
Creates an Upset plot to show overlap among the different samples, highlighting the peaks in gene bodies (using MaizeCode_R_Upset.r script)
Other analyses to add:
Calculates differential peaks between the different tissues (if the same mark is present in different tissue)
For RNA samples: NOT DONE YET, but expectations are:
Calls differentially expressed genes
Then, it combines both types of data: NOT DONE YET, but expectations are:
Plots heatmaps of all the samples over the regionfile (parameters will need to be worked on)
Plots metaplots of ChIP marks on the differentially expressed genes (each UP and DOWN between all sample pairs)
Plots Upset plots highlighting DEG in the overlapped peaks

- MaizeCode_R_Upset.r
Create an Upset plot of overlapping peaks and their presence in gene bodies (e.g. data/Upset_B73_endosperm.pdf)
