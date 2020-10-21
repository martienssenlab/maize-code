# maize-code

MaizeCode Pipeline Help


Step-by-Step pipeline

1) Make a MaizeCode folder somewhere
2) Check for packages requirements and install the ones you don’t have:
pigz; samtools; bowtie2; STAR; fastqc; cutadapt; bedtools; deeptools; macs2; idr;
3) Organize your reference genome folders so that they each contain ONE fasta file and ONE gff3 file
4) Copy all these scripts in ~/data/Scripts/. If you would rather have them in another folder, check in the MaizeCode.sh script and replace the location.
5) Make the samplefiles you want. Examples of how to make them are at the bottom of the MaizeCode.sh file (easiest is then to ctrl+H the names and path from the excel spreadsheet to replace them)
6) Submit the MaizeCode.sh script, giving as argument -f the samplefile.txt of your choice, -t the type of data (ChIP, RNA, ...), -p the path to your directory that contains the different genome directories and -r the name of the reference to use (i.e. the name of the genome directory folder).
7) Make the analysis_samplefiles you want. Examples are at the bottom of the MaizeCode_ChIP_analysis.sh file.
8) For ChIP data, submit the MaizeCode_ChIP_analysis.sh script either in the MaizeCode or in the ChIP folder (gets created in MaizeCode.sh), giving as argument -f the analysis_samplefile.txt and -r the regions (bed file) to be plotted on. -s can be set if the analysis should stop after calling peaks and making bigwig files.


Comments

- Check out the usage of each script before as it can give some information about little things that might need to be edited. Submitting each script without arguments (or followed by help) will return its usage.
- The RNA pipeline is not ready yet
- It should work for both Single-end or Paired-end data but the SE part has not been tested (and I might not have edited it well as I was changing the PE part). A non-issue for now since all the ChIP data is PE, but to keep in mind for potential future use.
- A lot of the steps are only happening in some cases that should be more or less self-explanatory.
- The whole pipeline creates a lot of report files and probably files that are not necessary to keep but for now I keep them like this for potential troubleshooting.
- For now I’ve used the MaizeCode.sh script for 8 samples at a time (all ChIPs from a specific LinexTissue). It runs in a couple hours (depending on the size of the files). It processes each sample in parallel and asks for 20 threads for each so better to limit the samplefiles to this amount.
- Always process the Input samples with their corresponding ChIP in the MaizeCode.sh script. (It can also be done separately, before or after) but they need to be done for the MaizeCode_ChIP_analysis.sh script to run successfully.
- The MaizeCode_ChIP_analysis.sh script is clearly not finished and will have to be tweaked for each analysis, but should give a first look at the data.
- These are preliminary version of the scripts, Evan will probably have some input on them, so they are not final versions! Use at your own risk (i.e. can be good for preliminary but probably not the final results). Comments from anyone are of course appreciated.
- I am still working on the MaizeCode_ChIP_analysis.sh script. The version here is far from the final version. It is just here for now to satisfy your curiosity. I’ll keep updating the file on the drive. The first part (calling peaks and making bigwig) should be good and you could run it if you want, but the rest of the analysis needs work. In any case, it will have to be tweaked for each analysis as we go, but I hope to make it good enough to run once for a first overview of all the data before specific analysis.


Scripts description

- MaizeCode.sh
Create the different folders
Create indexes (either bt2 or star) and reference files (if needed)
Copy fastq files from their folder in seq/ to the fastq/ folder (if not already done)
Run an instance of MaizeCode_ChIP_sample.sh or MaizeCode_RNA_sample.sh for each sample

-MaizeCode_ChIP_analysis.sh
Run fastQC on the raw data
Trim adapters, low quality and small reads (<20bp) with cutadapt
Run fastQC on trimmed data
Map with bowtie2
Removes PCR duplicates with samtools
Gets some mapping stats (in reports/summary_mapping_stats.txt)

-MaizeCode_ChIP_analysis.sh
Take either each replicate independently or merge the bam files of the two replicates
Call peaks with macs2 (calls broad peaks for H3K4me1, and narrow peaks for H3K4me3 and H3K27ac)
Make bigwig files with deeptools (log2 FC vs Input, normalizing each file by CPM) 
Do IDR analysis of the replicates with idr (if possible)
Get intersection of peaks to create different groups (presence/absence of each mark) with bedtools merge/intersect
Get distance of peaks in each group to closest gene with bedtools closest
Make heatmaps and metaplots with deeptools (on all genes (scaled by region) and around the TSS, and split into 5 clusters). This step will be very dependent on what bedfile you want to use, so it will probably need case-specific edits. 
