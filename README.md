# maize-code


## MaizeCode Pipeline Help


### Step-by-Step pipeline

1) Make a MaizeCode folder somewhere
2) Copy the MaizeCode scripts folder somewhere
2B) In MaizeCode.sh, MaizeCode_analysis.sh and MaizeCode_combined_analysis.sh, replace the path to the folder containing the MaizeCode scripts with yours:\
`export mc_dir="${HOME}/data/Scripts/MaizeCode/"` with\
`export mc_dir="${HOME}/YOUR/PATH/TO/SCRIPTS/MaizeCode/"` 
__This is one of the TODO things to improve!!__
3) Check that the following required packages are installed and in your $PATH (the versions noted here are working for sure, no guarantees for different versions). Recommended installation using conda.
```
pigz 2.3.4
samtools 1.10 (Using htslib 1.10.2)
bowtie2 64-bit 2.4.1; Compiler: gcc version 7.5.0 (crosstool-NG 1.24.0.131_87df0e6_dirty)
STAR 2.7.5c
fastqc v0.11.9
cutadapt 2.10
bedtools v2.29.2
deeptools 3.5.0
macs2 2.2.7.1
IDR 2.0.4.2 
R 3.6.3
R libraries: readr 1.4.0; ggplot2 3.3.2; UpSetR 1.4.0
```
4) Organize your reference genome folders so that they are all in the same main folder and that each contain ONE fasta file and ONE gff3 file (has to have 'gene' in column 3 and exons must be linked by 'Parent' in column 9)
5) Make the samplefiles you want. An example of a samplefile is in the data folder (B73_endosperm_samplefile.txt) and a quick way to make them is at the bottom of the `MaizeCode.sh` file. For cleaner naming purposes, use "\_samplefile.txt" as a suffix.
6) Submit the `MaizeCode.sh` script, giving as argument `-f <samplefile>` the samplefile.txt of your choice and `-p <path> ` the path to your directory that contains the different genome directories.
7) By default, it will proceed with the analysis. `-s` can be set so that it does not proceed with the analysis at all, and it will not do the combined analysis if several different references are being used for mapping
8) If the analysis has not proceeded or if you want to analyze different samples together, make the analysis_samplefile you want. An example of an analysis samplefile is in the data folder (B73_endosperm_analysis_samplefile.txt) and a quick way to make them is at the bottom of the `MaizeCode_analysis.sh` file. For better naming purposes, use "\_analysis_samplefile.txt" as a suffix.
9) Submit the `MaizeCode_analysis.sh` script, giving as argument `-f <analysisfile>` the analysis_samplefile.txt and `-r <bedfile>` the regions (bed file) to be plotted on. `-s` can be set if the combined analysis should not be performed (only calls peaks and makes bigwig files). It also stops there if the regionfile is missing.
10) Have a look at the results: mapping statistics, peak statistics and various plots (see Output below).

---

### Comments

- There is one wrapper script `MaizeCode.sh` that launches sub-scripts depending on what needs to be done.
- The `MaizeCode_analysis.sh` script is called by default by the `MaizeCode.sh` script but can be used seperately for any additional analysis.
- Check out the usage of each script before as it can give some information about the details and requirement for the scripts. Submitting each script without arguments (or followed by `-h`) will return its usage.
- The RNA pipeline is not ready yet
- It should work for both Single-end or Paired-end data but the SE part has not been tested (and I might not have edited it well as I was changing the PE part). A non-issue for now since all the ChIP data is PE, but to keep in mind for potential future use.
- The whole pipeline creates a lot of report files and probably files that are not necessary to keep but for now I keep them like this for potential troubleshooting.
- For now I’ve used the `MaizeCode.sh` script for 8 samples at a time (all ChIPs from a specific Line x Tissue). It runs in a couple hours (depending on the size of the files). It processes each sample in parallel and asks for 20 threads for each so better to limit the samplefiles to this amount.
- Always process the Input samples with their corresponding ChIP in the `MaizeCode.sh` script. (It can also be done separately, before or after) but they need to be done for the `MaizeCode_analysis.sh` script to run successfully.
- The `MaizeCode_combined_analysis.sh` script will have to be tweaked for each analysis, but a default use of it should give a first look at the data.
- These are preliminary version of the scripts!

---

### Scripts description

- __MaizeCode.sh__ - _wrapper script for the whole pipeline_\
Creates the different folders\
Starts the `MaizeCode_check_environment.sh` script for each environment (datatype * reference) that needs to be created\
Copies fastq files from their folder in seq/ to the fastq/ folder (if not already done)\
Runs an instance of `MaizeCode_ChIP_sample.sh` or `MaizeCode_RNA_sample.sh` for each sample\
Waits for the samples to be mapped\
Launches the `MaizeCode_analysis.sh` script if the `-s` argument (that stops after mapping) has not been given\
If mapping all the files to the same reference genome, it will use the all_genes.bed file created by the check_environment script\
Otherwise, it will not performed the combined analysis (it would not provide a regionfile)

- __MaizeCode_check_environment.sh__\
Checks if there is ONE fasta and ONE gff3 file in the reference folder (and unzip them if required)\
Makes a `chrom.sizes` file if not there (can be useful down the line for bedGraphtoBigWig for example)\
Makes a `all_genes.bed` file if not there (will be used for analysis/plots)\
Create the template for the stat files (in `[ChIP|RNA]/reports/summary_mapping_stats.txt`)\
Makes the bowtie2 or STAR indexes (for ChIP and RNA, respectively) if not already there

- __MaizeCode_ChIP_sample.sh__\
Runs fastQC on the raw data\
Trims adapters, low quality and small reads (<20bp) with cutadapt\
Runs fastQC on trimmed data\
Maps with bowtie2\
Removes PCR duplicates with samtools\
Gets some mapping stats (in `ChIP/reports/summary_mapping_stats.txt`)

- __MaizeCode_RNA_sample.sh__ - ___NOT DONE YET, but expectations are:___\
Runs fastQC on the raw data\
Trims adapters and low quality with cutadapt\
Runs fastQC on trimmed data\
Maps with STAR with different settings depending on the type of data (RNAseq, shRNA, RAMPAGE)\
Gets some mapping stats (in `RNA/reports/summary_mapping_stats.txt`)

- __MaizeCode_analysis.sh__ - _wrapper script for the analysis pipeline_\
If new files are to be analyzed, send each group of sample of the same datatype to `MaizeCode_ChIP_analysis.sh` or `MaizeCode_RNA_analysis.sh`\
Launches the `MaizeCode_combined_analysis.sh` if the `-s` option is not set (and if a regionfile is given)

- __MaizeCode_ChIP_analysis.sh__\
Merges biological replicates and split into pseudo-replicates\
For each type of file (replicate1, replicate2, pseudo-replicate1, pseudo-replicate2 and merged) in parallel:
  Calls peaks with macs2 (calls broad peaks for H3K4me1, and narrow peaks for H3K4me3 and H3K27ac)\
  Makes bigwig files with deeptools (log2 FC vs Input, normalizing each file by CPM)\
  Plot Fingerprint\
Waits for the previous steps to proceed
Makes IDR analysis for biological replicates with idr\
Makes a `selected_peaks` file with the peaks called in the merged sample and both pseudo-replicates with bedtools intersect\
Makes some stats on the number of peaks (in `ChIP/peaks/summary_peaks_<samplefile_name>.txt`)

- __MaizeCode_RNA_analysis.sh__ - ___NOT DONE YET, but expectations are:___\
Makes bigwig files with deeptools (normalized by CPM)
Calls expressed genes for each sample (_min CPM to be determined_)

- __MaizeCode_combined_analysis.sh__\
Splits the samplefile into ChIP and RNA samples\
For ChIP samples:\
Makes a single file, merging all selected peaks from all samples with bedtools merge\
Gets distance of each peak to the closest region from the regionfile with bedtools closest (default: all genes annotated in the reference)\
Creates an Upset plot to show overlap among the different samples, highlighting the peaks in gene bodies (using `MaizeCode_R_Upset.r` script)\
___Other analyses to add:___\
Calculates differential peaks between the different tissues (if the same mark is present in different tissue) _in progress_\
For RNA samples: ___NOT DONE YET, but expectations are:___\
Calls differentially expressed genes\
Then, it combines both types of data: ___NOT DONE YET, but expectations are:___\
Plots heatmaps of all the samples over the regionfile (parameters will need to be worked on)\
Plots metaplots of ChIP marks on the differentially expressed genes (each UP and DOWN between all sample pairs)\
Plots Upset plots highlighting expressed genes in the overlapped peaks

- __MaizeCode_R_Upset.r__\
Create an Upset plot of overlapping peaks and their presence in gene bodies. An example is in the data folder (e.g. data/Upset_B73_endosperm.pdf)
___Other option to add:___\
Highlight expressed genes instead of gene bodies

---

### Output

__Directories:__
From the main folder `<maizecode>` where the `MaizeCode.sh` is run

- `<maizecode>/ChIP`: Folder containing data from ChIP sample(s)
*only created if at least one ChIP sample has been analyzed*
  - `<maizecode>/ChIP/fastq`: Folder containing raw and trimmed fastq files
  - `<maizecode>/ChIP/mapped`: Folder containing mapped and indexed data (bam and bam.bai files). It will contain mapped data before and after deduplication for each biological replicate, the merged replicates and the pseudo-replicates files.
  - `<maizecode>/ChIP/tracks`: Folder containing bigwig files and the all_genes.bed file for all genome references
  - `<maizecode>/ChIP/plots`: Folder containing the fingerprint plots for each sample and idr plots between biological replicates
  - `<maizecode>/ChIP/peaks`: Folder containing all peak files and the `summary_peaks_<samplefile_name>.txt` file that has a summary of peak statistics for all ChIP samples analyzed together in the samplefile `<samplefile_name>_samplefile.txt`
  - `<maizecode>/ChIP/reports`: Folder containing the fastQC reports, trimming details, mapping details, idr details and the `summary_mapping_stats.txt` file that has a summary of mapping statistics for all ChIP samples processed
  - `<maizecode>/ChIP/logs`: Folder containing log files to go back to in case of error during environment building `env_<genome_reference>.log`, mapping `<sample_name>.log`, analysis of ChIP samples together `<samplefile_name>.log` and single sample analysis `analysis_<sample_name>.log`
  - `<maizecode>/ChIP/chkpts`: Folder containing `touch` files to track success and completion of environment building `env_<genome_ref>`, sample mapping `<sample_name>` and single sample analysis (peak calling and bigwig files) `analysis_<sample_name>`. These files are produced to prevent these steps to be repeated if they were already performed in order to only performed the combined analysis of different combinations of samples. If these files are deleted, the mapping and analysis steps will be repeated and will overwrite existing files.

- `<maizecode>/RNA`: Folder containing data from RNA sample(s) ___NOT DONE YET, but expectations are:___\
*only created if at least one ChIP sample has been analyzed*
  - `<maizecode>/RNA/fastq`: Folder containing raw and trimmed fastq files
  - `<maizecode>/RNA/mapped`: Folder containing mapped and indexed data (bam and bam.bai files). It will contain mapped data before and after deduplication for each biological replicate and the merged replicates files.
  - `<maizecode>/RNA/tracks`: Folder containing bigwig files and the all_genes.bed file for all genome references
  - `<maizecode>/RNA/plots`: Folder containing some plots for each sample (*which plots, if any, to be determined*)
   - `<maizecode>/RNA/reports`: Folder containing the fastQC reports, trimming details, mapping details and the `summary_mapping_stats.txt` file that has a summary of mapping statistics for all RNA samples processed
  - `<maizecode>/RNA/logs`: Folder containing log files to go back to in case of error during environment building `env_<genome_reference>.log`, mapping `<sample_name>.log` and analysis of RNA samples together `<samplefile_name>.log`
  - `<maizecode>/RNA/chkpts`: Folder containing `touch` files to track success and completion of environment building `env_<genome_ref>`, sample mapping `<sample_name>` and single sample analysis (bigwig files) `analysis_<sample_name>`. These files are produced to prevent these steps to be repeated if they were already performed in order to only performed the combined analysis of different combinations of samples. If these files are deleted, the mapping and analysis steps will be repeated and will overwrite existing files.

- `<maizecode>/combined`: Folder containing data from combined analysis ___The names of the folders should be changed to be more explicit___\
*only created if at least one combined analysis has been performed*
  - `<maizecode>/combined/DEG`: Folder containing differentially expressed genes analysis results
  - `<maizecode>/combined/peaks`: Folder containing combined peak files and matrix for Upset plots
  - `<maizecode>/combined/matrix`: Folder containing matrix files for heatmap plotting, outputed regions from kmean clustering of the heatmaps and value tables to be used for the scales of heatmaps
  - `<maizecode>/combined/plots`: Folder containing the Upset plots and heatmaps
  - `<maizecode>/combined/logs`: Folder containing log files to go back to in case of error during combined analysis `combined_analysis_<samplefile_name>_<regionfile_name>.log`
  - `<maizecode>/combined/chkpts`: Folder containing `touch` files to track success of combined analysis `<samplefile_name>_<regionfile_name>.log`. These files are only for success tracking and will be overwritten if an analysis with the same name is to be performed. 

- `<maizecode>/chkpts`: Folder containing `touch` files to track success of a run **without** combined analysis

__Statistics:__
- `summary_mapping_stats.txt`
Located in `<maizecode>/ChIP/reports/` for ChIP samples and `<maizecode>/RNA/reports` for RNA samples.
Tab-delimited file with 8 columns giving information for each sample (detailed in columns#1 to #4) on\
the number of reads in total (column #5),\
the number of reads (and percentage of the total reads) that pass filtering (column #6),\
the number of reads (and percentage of the total reads) that are in the deduplicated bam file (column #7)\,
the number of reads (and percentage of the total reads) that are properly mapped (column #8).

- `summary_peaks_<samplefile_name>.txt`
Located in `<maizecode>/ChIP/peaks/`. 
Tab-delimited file with 10 columns giving information for each histone mark (detailed in columns#1 to #3) on\
the number of peaks called in each biological replicate (columns #4 and #5, respectively),\
the number of peaks in common between the biological replicates (all peaks given by the IDR analysis) and the percentage relative to each biological replicate (column #6),\
the number of peaks in common between the biological replicates that pass the IDR threshold of 0.05 and the percentage relative to the number of peaks in common (column #7),\
the number of peaks called when both replicates are merged (column #8),\
the number of peaks shared by each pseudo-replicate (column #9),\
the number of selected peaks (i.e. the peaks that will be used for downstream analysis) which are the peaks shared by the merged and both pseudo-replicates, and the percentage relative to the the number of merged peaks (column #10).

__Plots:__
- `Upset_<samplefile_name>_<regionfile_name>.pdf`
Upset plots showing intersection between all the samples in the `<samplefile_name>_analysis_samplefile.txt`, highlighting the peaks that are present on the regions in `<regionfile_name>.bed`

- `<samplefile_name>_<regionfile_name>_heatmap_regions.pdf`
Heatmap of the enrichment for all ChIP and RNA samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, scaling each region to the same length, in decreasing order of overall enrichment in all samples

- `<samplefile_name>_<regionfile_name>_heatmap_regions_k5.pdf`
Heatmap of the enrichment for all ChIP and RNA samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, scaling each region to the same length, clustered into 5 regions by kmeans

- `<samplefile_name>_<regionfile_name>_heatmap_tss.pdf`
Heatmap of the enrichment for all ChIP and RNA samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, aligning all regions by their transcription start site, in decreasing order of overall enrichment in all samples

- `<samplefile_name>_<regionfile_name>_heatmap_tss_k5.pdf`
Heatmap of the enrichment for all ChIP and RNA samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, aligning all regions by their transcription start site, clustered into 5 regions by kmeans

- `<samplefile_name>_<regionfile_name>_heatmap_tss.pdf`
Heatmap of the enrichment for all ChIP and RNA samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, aligning all regions by their transcription end site, in decreasing order of overall enrichment in all samples

- `<samplefile_name>_<regionfile_name>_heatmap_tes_k5.pdf`
Heatmap of the enrichment for all ChIP and RNA samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, aligning all regions by their transcription end site, clustered into 5 regions by kmeans

