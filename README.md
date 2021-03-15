# maize-code


## MaizeCode Pipeline Help


### Step-by-Step pipeline


1) Clone the git repository anywhere you want, for example in a new folder called projects\
`git clone https://github.com/eernst/maize-code.git ./projects/maize-code`\
or to clone a specific branch 'devel'\
`git clone --branch devel https://github.com/eernst/maize-code.git ./projects/maize-code`\
You will be prompted to input your GitHub username and password.\
If you only want to update the scripts, use `git pull`. If you want to update from a specific branch 'devel' `git pull origin devel`.
2) cd into the maize-code folder that has been created, so following the same example\
`cd ./projects/maize-code/`
3) Check that the following required packages are installed and in your $PATH (the versions noted here are working for sure, no guarantees for different versions). Recommended installation using conda (except grit that should be installed with pip, but finding an alternative to using it is being looked at)
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
(grit 2.0.5)
bedGraphToBigWig v 2.8
bedSort
parallel-fastq-dump (if downloading from SRA)
R 3.6.3
R libraries: ggplot2 3.3.2; UpSetR 1.4.0; limma 3.42.2; edgeR 3.28.1; dplyr 1.0.2; tidyr 1.1.2; stringr 1.4.0; cowplot 1.1.0; gplots 3.1.0; RColorBrewer 1.1.2; ComplexUpset ; purrr ; AnnotationForge ; rrvgo ;
```
4) Organize your reference genome directories so that they are all in the same main folder and that each contain ONE fasta file (.fa extension), ONE GFF file (.gff or .gff* extension) and ONE GTF (.gtf extension) file.\ 
For example, having a `genomes/` folder that contains the `genomes/B73/` directory where you can find `genomes/B73/B73.fa`, `genomes/B73/B73.gff` and `genomes/B73/B73.gtf` files\
Other references should be in the same `genomes/` folder, following the same pattern, i.e. `genomes/W22/W22.fa`, `genomes/W22/W22.gff` and `genomes/W22/W22.gtf`\
The GTF file can be created from a GFF file with cufflinks `gffread -T <gff_file> -o <gtf_file>` and check that 'transcript_id' and 'gene_id' look good in the 9th column.\
The GFF file should have 'gene' in the 3rd column.\
All files can be gzipped (.gz extension).
5) Make the samplefiles you want. An example of a samplefile is in the data folder (Example_samplefile.txt) and a quick way to make them is at the bottom of the `MaizeCode.sh` script. For cleaner naming purposes, use "\_samplefile.txt" as a suffix.
6) Submit the `scripts/MaizeCode.sh` script, giving as argument `-f <samplefile.txt>` the samplefile of your choice and `-p <path>` the path to your folder that contains the different genome directories, i.e. the `genomes` folder mentioned above:\
`qsub scripts/MaizeCode.sh -f example_samplefile.txt -r /path/to/genomes`
7) By default, it will proceed with the analysis. `-s` can be set so that it does not proceed with the analysis at all or `-c` can be set if only single sample analysis should be performed but no combined analysis per line or between lines.
8) If the analysis has not proceeded or if you want to analyze different samples together, make the new samplefile of your choice and submit the `scripts/MaizeCode.sh` script again.\
`qsub scripts/MaizeCode.sh -f new_samplefile.txt -r /path/to/main/folder/containing/genome/directories`\
The samples that have already been processed will not be repeated but will still be included in the analysis.
9) Have a look at the results! (see Output below).

---

### Comments

- There is one wrapper script `MaizeCode.sh` that launches sub-scripts depending on what needs to be done.
- Potentially, each script could be submitted on its own but it could be tricky. Check out the usage of each script before by running the script without arguments (or followed by `-h`).
- The shRNA and complete RAMPAGE pipelines are not ready yet
- It should work for both Single-end or Paired-end data
- It works perfectly with 2 replicates for every type of data (including two different inputs for ChIP). Adapting the scripts to allow for more variation in the number of replicates is under development (having only one ChIP input replicate works, as well as multiple RNAseq replicates).
- The whole pipeline creates a lot of report files and probably files that are not necessary to keep but for now I keep them like this for potential troubleshooting.
- For now I’ve used the `MaizeCode.sh` script from scratch for all samples of one inbred line at a time (all ChIPs and RNAseq from 5 tissues). It runs in ~19h (but it depends on the size of the files and how busy the cluster is...). Once that the mapping and single-sample analysis have been done, reusing these samples in a different analysis is much quicker though (~2h), the limitations are for mapping ChIPseq samples and calling ChIPseq peaks (since it does it for each biological replicate, the merge file and both pseudo-replicates and cannot be multi-threaded). That is probably the first step that could be optimized for faster runs.
- Always process the Input samples with their corresponding ChIP in the `MaizeCode.sh` script.
- The analysis will have to be adapted to the desired output, but running the default complete pipeline should give a first look at the data and generate all the files required for further analysis.
- These are still preliminary version of the scripts!

---

### Scripts description

- __MaizeCode.sh__ - _wrapper script for the whole pipeline_\
Creates the different folders\
Runs the `MaizeCode_check_environment.sh` script for each environment (datatype * reference) that needs to be created\
Runs an instance of `MaizeCode_ChIP_sample.sh` or `MaizeCode_RNA_sample.sh` for each sample\
Waits for the samples to be mapped\
Runs the `MaizeCode_R_mapping_stats.r` script to plot the mapping statitistics of all the samples in samplefile into bar plots\ 
Runs the `MaizeCode_analysis.sh` script if the `-s` argument (that stops after mapping) has not been given\
Runs the `MaizeCode_analysis.sh -s` script if the `-c` argument (that stops after single sample analyis) has been given\
By default, it will provide the `<reference_genome>_all_genes.bed` files created by the check_environment script as region files

- __MaizeCode_check_environment.sh__\
Checks if there is ONE fasta and ONE gff3 file in the reference folder (and unzip them if required)\
Makes a `chrom.sizes` file if not there (can be useful down the line for bedGraphtoBigWig for example)\
Makes a `all_genes.bed` file if not there (will be used for analysis/plots)\
Create the template for the stat files\
Makes the bowtie2 or STAR indexes (for ChIP and RNA, respectively) if not already there

- __MaizeCode_ChIP_sample.sh__\
Copies fastq files from their original folder or GEO to the fastq/ folder (if not already done)\
Runs fastQC on the raw data\
Trims adapters, low quality and small reads (<20bp) with cutadapt\
Runs fastQC on trimmed data\
Maps with bowtie2\
Removes PCR duplicates with samtools\
Gets some mapping stats

- __MaizeCode_RNA_sample.sh__ ___shRNA NOT DONE YET and will be processed by a different pipeline (using shortStack)___\
Copies fastq files from their original folder or GEO to the fastq/ folder (if not already done)\
Runs fastQC on the raw data\
Trims adapters and low quality with cutadapt\
Runs fastQC on trimmed data\
Maps with STAR with different settings depending on the type of data (RNAseq, shRNA, RAMPAGE)\
Marks duplicates with STAR with different settings depending on the type of data (RNAseq, shRNA, RAMPAGE)\
Creates stranded bigwig files with STAR and bedGraphToBigWig with different settings depending on the type of data (RNAseq, shRNA, RAMPAGE)\
Gets some mapping stats

- __MaizeCode_analysis.sh__ - _wrapper script for the analysis pipeline_\
If new samples are to be analyzed individually, sends each group of samples of the same datatype to `MaizeCode_ChIP_analysis.sh` or `MaizeCode_RNA_analysis.sh`\
Gathers peak statistics for all ChIPseq samples in the samplefile and lauches `MaizeCode_R_peak_stats.r` to plot them\
Gathers gene expression statistics for all RNAseq samples in the samplefile and lauches `MaizeCode_R_gene_ex_stats.r` to plot them\
Gathers TSS statistics for all RAMPAGE samples in the samplefile\
Launches the `MaizeCode_line_analysis.sh` script for each reference present in the samplefile\
If different lines are present, launches the `MaizeCode_combined_analysis.sh` script ___NOT DONE YET___

- __MaizeCode_ChIP_analysis.sh__\
Merges biological replicates and split into pseudo-replicates\
For each type of file (replicate1, replicate2, pseudo-replicate1, pseudo-replicate2 and merged) in parallel:\
Calls peaks with macs2 (calls broad peaks for H3K4me1, and narrow peaks for H3K4me3 and H3K27ac)\
Makes bigwig files with deeptools (log2 FC vs Input, normalizing each file by CPM)\
Plot Fingerprint\
Waits for the previous steps to proceed\
Makes IDR analysis for biological replicates with idr\
Makes a `selected_peaks` file with the peaks called in the merged sample and both pseudo-replicates with bedtools intersect\
Makes some stats on the number of peaks

- __MaizeCode_RNA_analysis.sh__\
Processes each sample in parallel\
For RAMPAGE data:\
Merges biological replicates and creates stranded tracks (bigwigs) with STAR and bedGraphToBigWig\
Calls peaks (to identify TSS) with macs2 (_should be grit but not maintained and pretty cryptic_)\
Makes IDR analysis for biological replicates with idr\
Make some stats on the number of peaks/tss\
For RNAseq data:\
Merges biological replicates and creates stranded tracks (bigwigs) with STAR and bedGraphToBigWig\
Makes some stats on the number of expressed genes\

- __MaizeCode_line_analysis.sh__ ___Analyses marked by *** are still under development:___\
Splits the samplefile into ChIPseq and RNA samples\
For ChIPseq samples:\
Makes a single file, merging all selected peaks from all samples with bedtools merge\
Gets distance of each peak to the closest region from the regionfile with bedtools closest (default: all genes annotated in the reference)\
Creates an Upset plot to show overlap among the different samples, highlighting the peaks in gene bodies, using `MaizeCode_R_Upset.r` script\
_if several tissues are present in the samplefile:_\
Calculates differential peaks between the different tissues\*\*\*\
For RAMPAGE samples:\
_if several tissues are present in the samplefile:_\
Calls differential TSS between the different tissues\*\*\*\
For RNAseq samples:\
_if several tissues are present in the samplefile:_\
Makes sample and count tables\
Calls differentially expressed genes between all pairs of tissues using `MaizeCode_R_DEG.r` script\
On all the samples:
Plots heatmaps of the ChIPseq and RNAseq samples over the regionfile (parameters might need to be adjusted)\
Plots heatmaps and profiles of the ChIPseq samples over the differentially expressed genes (if they were called previously)\
Splits all genes into different clusters based on all avalaible data (silent, constitutive and tissue-specific genes)\*\*\*\
Identifies enhancers and assign to a gene\*\*\*

- __MaizeCode_combined_analysis.sh__ ___NOT DONE YET, but expectations are:___\
Compares gene status between homolog genes\
Compares enhancers?

- __MaizeCode_R_mapping_stats.r__\
Creates a plot representing the mapping statistics (both read numbers and distribution), named `combined/plots/mapping_stats_<samplefile_name>.pdf`

- __MaizeCode_R_peak_stats.r__\
Creates a plot representing the number of peaks called in the different ChIPseq sub-samples, named `combined/plots/peak_stats_<samplefile_name>.pdf`

- __MaizeCode_R_gene_ex_stats.r__\
Creates a plot summarizing the expression levels of genes in all RNAseq samples named `combined/plots/gene_expression_stats_<samplefile_name>.pdf`

- __MaizeCode_R_Upset.r__\
Creates an Upset plot of overlapping peaks and their presence in gene bodies named `combined/plots/Upset_<analysis_name>.pdf`

- __MaizeCode_R_DEG.r__\
Performs differential expression analysis with edgeR on all RNAseq samples present in the samplefile with edgeR\
Plots MDS (`combined/plots/MDS_<analysis_name>.pdf`) and BCV (`combined/plots/BCV_<analysis_name>.pdf`)\
For each pair of tissues:\
Create a table of log2FC for all genes (named `combined/DEG/FC_<analysis_name>_<tissue1>_vs_<tissue2>.txt`)\
Create a table of differentially expressed genes (named `combined/DEG/DEG_<analysis_name>_<tissue1>_vs_<tissue2>.txt`)\
Plots two heatmaps on all the differentially expressed genes (by log(cpm) named `combined/plots/Heatmap_cpm_<analysis_name>.pdf` and scaling per row (z_score) named `combined/plots/Heatmap_zscore_<analysis_name>.pdf`)

---

### Output

__Directories:__
From the main folder `<maizecode>` where the `MaizeCode.sh` is run

- `<maizecode>/ChIP`: Folder containing data from ChIPseq sample(s)\
*only created if at least one ChIPseq sample has been analyzed*
  - `<maizecode>/ChIP/fastq`: Folder containing raw and trimmed fastq files
  - `<maizecode>/ChIP/mapped`: Folder containing mapped and indexed data (bam and bam.bai files). It will contain mapped data before and after deduplication for each biological replicate, the merged replicates and the pseudo-replicates files.
  - `<maizecode>/ChIP/tracks`: Folder containing bigwig files and the all_genes.bed file for all genome references
  - `<maizecode>/ChIP/plots`: Folder containing the fingerprint plots for each sample and idr plots between biological replicates
  - `<maizecode>/ChIP/peaks`: Folder containing all peak files and idr analysis between biological replicates
  - `<maizecode>/ChIP/reports`: Folder containing the fastQC reports, trimming details, mapping details, idr details, the `summary_mapping_stats.txt` file that has a summary of mapping statistics for all ChIPseq samples ever processed and the `summary_ChIP_peaks.txt` file that has a summary of peak statistics for all ChIPseq samples ever processed
  - `<maizecode>/ChIP/logs`: Folder containing log files to go back to in case of error during environment building `env_<genome_reference>.log`, mapping `<sample_name>.log`, analysis of ChIP samples together `<samplefile_name>.log` and single sample analysis `analysis_<sample_name>[|_Rep1|_Rep2|_pseudo1|_pseudo2|_merged].log`
  - `<maizecode>/ChIP/chkpts`: Folder containing `touch` files to track success and completion of environment building `env_<genome_ref>`, sample mapping `<sample_name>` and single sample analysis (peak calling and bigwig files) `analysis_<sample_name>`. These files are produced to prevent these steps to be repeated if they were already performed in order to only performed the combined analysis of different combinations of samples. If these files are deleted, the mapping and analysis steps will be repeated and will overwrite existing files.

- `<maizecode>/RNA`: Folder containing data from RNA sample(s)\
*only created if at least one RNA sample has been analyzed*
  - `<maizecode>/RNA/fastq`: Folder containing raw and trimmed fastq files
  - `<maizecode>/RNA/mapped`: Folder containing mapped and indexed data (bam and bam.bai files). It will contain mapped data before and after deduplication for each biological replicate and the merged replicates files
  - `<maizecode>/RNA/tracks`: Folder containing stranded bigwig files based on all or unique reads (4 files per biological replicate, plus 4 files for the merged replicates) and the all_genes.bed file for all genome references
  - '<maizecode>/RNA/TSS`: Folder containing the peaks (TSS) called on RAMPAGE data and the IDR analysis between biological replicates
  - `<maizecode>/RNA/plots`: Folder containing the idr plots between biological replicates for RAMPAGE samples
  - `<maizecode>/RNA/reports`: Folder containing the fastQC reports, trimming details, mapping details, the `summary_mapping_stats.txt` file that has a summary of mapping statistics for all RNA samples ever processed, the `summary_gene_expression.txt` file that has a summary of gene expression statistics for all RNAseq samples ever processed and the `summary_RAMPAGE_tss.txt` file that has a summary of peak/tss statistics for all RAMPAGE samples ever processed
  - `<maizecode>/RNA/logs`: Folder containing log files to go back to in case of error during environment building `env_<genome_reference>.log`, mapping `<sample_name>.log`, analysis of RNA samples together `<samplefile_name>.log` and single sample analysis `analysis_<sample_name>.log`
  - `<maizecode>/RNA/chkpts`: Folder containing `touch` files to track success and completion of environment building `env_<genome_ref>`, sample mapping `<sample_name>` and single sample analysis (bigwig files) `analysis_<sample_name>`. These files are produced to prevent these steps to be repeated if they were already performed in order to only performed the combined analysis of different combinations of samples. If these files are deleted, the mapping and analysis steps will be repeated and will overwrite existing files.

- `<maizecode>/combined`: Folder containing data from combined analysis\
*only created if at least one sample has been mapped*\
**`<analysis_name>` is a combination of the samplefile and regionfile names: `<samplefile_name>_on_<regionfile_name>`**\
**By default, `<regionfile_name>` is `all_genes`**
  - `<maizecode>/combined/DEG`: Folder containing differentially expressed genes analysis results. `FC_<sample1>_vs_<sample2>.txt` are pairwise comparison between sample1 and sample2 for all genes. `DEG_<sample1>_vs_<sample2>.txt` only contain the differentially expressed genes (FDR<=0.05, |logFC|>2) between sample1 and sample2.
  - `<maizecode>/combined/peaks`: Folder containing combined ChIPseq peak files `peaks_<analysis_name>.bed`, RAMPAGE TSS files `tss_<analysis_name>.bed` and matrix for Upset plots `matrix_upset_<analysis_name>.txt`.
  - `<maizecode>/combined/matrix`: Folder containing matrix files for heatmap plotting `regions_<analysis_name>.gz`, `tss_<analysis_name>.gz` and `deg_<analysis_name>.gz`, outputed regions from kmean clustering of the heatmaps `<analysis_name>_regions_regions_k5.txt` and `<analysis_name>_tss_regions_k5.txt` and the value tables to be used for the scales of heatmaps `values_regions_<analysis_name>.txt` and `values_tss_<analysis_name>.txt`. 'regions' corresponds to the 'scale_regions' argument of deeptools, 'tss' corresponds to the 'reference-point --referencePoint TSS' argument of deeptools and 'k5' corresponds to the '--kmeans 5' argument of deeptools.
  - `<maizecode>/combined/plots`: Folder containing the mapping statistics plot `mapping_stats_<analysis_name>.pdf`, the peak statistics plot `peak_stats_<analysis_name>.pdf`, the gene expression statistics plot `gene_expression_stats_<analysis_name>.pdf`, the Upset plots of peaks in gene bodies `Upset_<analysis_name>.pdf`, the MDS and BCV plots from the DEG analysis `MDS_<analysis_name>.pdf` and `BCV_<analysis_name>.pdf`, respectively, the heatmaps of differentially expressed genes clustered accross all samples with log(cpm) values `Heatmap_cpm_<analysis_name>.pdf` and normalized for each gene `Heatmap_zscore_<analysis_name>.pdf`, and the different deeptools heatmaps `<analysis_name>_heatmaps_regions.pdf`, `<analysis_name>_heatmaps_regions_k5.pdf`, `<analysis_name>_heatmaps_tss.pdf` and `<analysis_name>_heatmaps_tss_k5.pdf`, and profiles . 'regions' corresponds to the 'scale_regions' argument of deeptools, 'tss' corresponds to the 'reference-point --referencePoint TSS' argument of deeptools and 'k5' corresponds to the '--kmeans 5' argument of deeptools.
  - `<maizecode>/combined/reports`: Folder containing the `summary_mapping_stats_<samplefile_name>.txt` file that has a summary of mapping statistics for all samples in each samplefile, the `summary_ChIP_peaks_<samplefile_name>.txt` file that has a summary of peak statistics for all ChIPseq samples in each samplefile,the `summary_gene_expression_<samplefile_name>.txt` file that has a summary of gene expression statistics for all RNAseq samples in each samplefile and the `summary_RAMPAGE_tss_<samplefile_name>.txt` file that has a summary of peak/tss statistics for all RAMPAGE samples in each samplefile.
  - `<maizecode>/combined/logs`: Folder containing log files to go back to in case of error during combined analysis `analysis_<analysis_name>.log`
  - `<maizecode>/combined/chkpts`: Folder containing `touch` files to track success of combined analysis `<analysis_name>.log`. These files are only for success tracking and will be overwritten if an analysis with the same name is to be performed.

- `<maizecode>/chkpts`: Folder containing `touch` files to track success of a run **without** analysis

__Statistics:__
- `summary_mapping_stats.txt`
Located in `<maizecode>/<ChIP|RNA>/reports/`\
Tab-delimited file with 8 columns giving information for each sample (detailed in columns#1 to #4) on\
the genome reference it was mapped to (column #5),\
the number of reads in total (column #6),\
the number of reads (and percentage of the total reads) that pass filtering (column #7),\
the number of reads (and percentage of the total reads) that are mapping to the reference (inlcuding multi-mappers) (column #8)\,
the number of reads (and percentage of the total reads) that are uniquely mapped (column #9).

- `summary_ChIP_peaks.txt` and `summary_ChIP_peaks_<samplefile_name>.txt`
Located in `<maizecode>/combined/reports/` and `<maizecode>/ChIP/reports/`\
Tab-delimited file with 10 columns giving information for each histone mark (detailed in columns#1 to #3) on\
the number of peaks called in each biological replicate (columns #4 and #5, respectively),\
the number of peaks in common between the biological replicates (all peaks given by the IDR analysis) and the percentage relative to each biological replicate (column #6),\
the number of peaks in common between the biological replicates that pass the IDR threshold of 0.05 and the percentage relative to the number of peaks in common (column #7),\
the number of peaks called when both replicates are merged (column #8),\
the number of peaks shared by each pseudo-replicate (column #9),\
the number of selected peaks (i.e. the peaks that will be used for downstream analysis) which are the peaks shared by the merged and both pseudo-replicates, and the percentage relative to the the number of merged peaks (column #10).

- `summary_RAMPAGE_tss.txt` and `summary_RAMPAGE_tss_<samplefile_name>.txt`
Located in `<maizecode>/combined/reports/` and `<maizecode>/RNA/reports/`\
Tab-delimited file with 8 columns giving information for each RAMPAGE sample (detailed in columns#1 to #3) on\
the number of annotated genes in the reference genome (columns #4),\
the number of peaks called in each biological replicate (columns #5 and #6, respectively),\
the number of peaks in common between the biological replicates (all peaks given by the IDR analysis) and the percentage relative to each biological replicate (column #7),\
the number of peaks in common between the biological replicates that pass the IDR threshold of 0.05 and the percentage relative to the number of peaks in common (column #8).\

- `summary_gene_expression.txt` and `summary_gene_expression_<samplefile_name>.txt`
Located in `<maizecode>/combined/reports/` and `<maizecode>/RNA/reports/`\
Tab-delimited file with 13 columns giving information for each RNAseq sample (detailed in columns#1 to #3) on\
the number of annotated genes in the reference genome (columns #4),\
the number of silent genes (cpm=0), lowly expressed genes (cpm<1) and highly expressed genes (cpm>1) in the first biological replicate (columns #5, #6 and #7, respectively),\ 
the number of silent genes (cpm=0), lowly expressed genes (cpm<1) and highly expressed genes (cpm>1) in the second biological replicate (columns #8, #9 and #10, respectively)\
the number of silent genes (cpm=0), lowly expressed genes (cpm<1) and highly expressed genes (cpm>1) after averaging both replicates (columns #11, #12 and #13, respectively)\

__Plots:__ (examples are in the github data folder)
- `ChIP/plots/Fingerprint_<sample_name>_<replicate>.png`\
Fingerprint plot from deeptools to assess the genome-wide distribution of reads for each ChIPseq sample replicate (Rep1, Rep2, merged, pseudo1 and pseudo2) and its corresponding Input.\
Tool details: https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html \
\[ generated by deeptools from `MaizeCode_ChIP_analysis.sh` ]

- `ChIP/plots/idr_<ChIPseq_sample_name>.png`\
Scatter plots and box plots from idr showing correlation between biological replicates of a ChIPseq sample.\
Tool details: https://github.com/nboley/idr \
\[ generated by idr from `MaizeCode_ChIP_analysis.sh` ]

- `RNA/plots/idr_<RAMPAGE_sample_name>.png`\
Scatter plots and box plots from idr showing correlation between biological replicates of a RAMPAGE sample.\
Tool details: https://github.com/nboley/idr \
\[ generated by idr from `MaizeCode_RNA_analysis.sh` ]

- `combined/plots/mapping_stats_<samplefile_name>.pdf`\
Bar plots showing number and distribution of uniquely mapped, multi-mapping, unmapped and filtered reads for all samples in `<samplefile_name>_analysis_samplefile.txt`.\
tool details: https://ggplot2.tidyverse.org/reference/geom_bar.html \
tool details: https://wilkelab.org/cowplot/reference/index.html \
\[ generated in R by `MaizeCode_R_mapping_stats.r`, started from `MaizeCode.sh` ]

- `combined/plots/peak_stats_<samplefile_name>.pdf`\
Bar plots showing the number of peaks called in each biological replicate, the number of peaks in common between the replicates, the number of peaks in common passing an IDR threshold of 0.05, the number of peaks called when the bam files of both biological replicates are merged, the number of peaks in common between pseudo-replicates (called on two random halves of the merged bam file) and the number of selected peaks (common between the merged and the pseudoreplicates) for all ChIPseq samples in `<samplefile_name>_analysis_samplefile.txt`.\
tool details: https://ggplot2.tidyverse.org/reference/geom_bar.html \
\[ generated in R by `MaizeCode_R_mapping_stats.r`, started from `MaizeCode.sh` ]

- `combined/plots/gene_expression_stats_<samplefile_name>.pdf`\
Bar plots showing the distribution of genes in three large categories: unexpressed (cpm=0), lowly expressed (cpm<1) and highly expressed (cpm>1) for each biological replicate and when taking their average, for all RNAseq samples in `<samplefile_name>_analysis_samplefile.txt`.\
tool details: https://ggplot2.tidyverse.org/reference/geom_bar.html \
\[ generated in R by `MaizeCode_R_mapping_stats.r`, started from `MaizeCode.sh` ]

- `combined/plots/Upset_<samplefile_name>_on_<regionfile_name>.pdf`\
Upset plots showing intersection between all the ChIP samples in the `<samplefile_name>_analysis_samplefile.txt`, highlighting the peaks that are present on the regions in `<regionfile_name>.bed`.\
Tool details: https://github.com/hms-dbmi/UpSetR \
\[ generated in R by `MaizeCode_R_Upset.r`, started from `MaizeCode_line_analysis.sh` ]

- `combined/plots/MDS_<samplefile_name>_on_<regionfile_name>.pdf`\
MDS plot (2D representation of variance) between all the RNAseq samples in the `<samplefile_name>_analysis_samplefile.txt` mapping to the same reference genome used to create the `<regionfile_name>.txt`.\
Tool details: https://rdrr.io/bioc/edgeR/man/plotMDS.DGEList.html \
\[ generated in R by `MaizeCode_R_DEG.r`, started from `MaizeCode_line_analysis.sh` ]

- `combined/plots/BCV_<samplefile_name>_on_<regionfile_name>.pdf`\
BCV plot (Biological coefficient of variation) for all the genes in the `<regionfile_name>.txt` based on all the RNAseq samples in the `<samplefile_name>_analysis_samplefile.txt`.\
Tool details: https://rdrr.io/bioc/edgeR/man/plotBCV.html \
\[ generated in R by `MaizeCode_R_DEG.r`, started from `MaizeCode_line_analysis.sh` ]

- `combined/plots/Heatmap_cpm_<samplefile_name>_on_<regionfile_name>.pdf`\
Clustered heatmap of all the differentially expressed genes between all pairs of RNAseq samples in the `<samplefile_name>_analysis_samplefile.txt`, scaled by log(count per million) of the RNAseq replicate samples (highlights gene expression levels).\
Tool details: https://www.rdocumentation.org/packages/gplots/versions/3.1.0/topics/heatmap.2 \
\[ generated in R by `MaizeCode_R_DEG.r`, started from `MaizeCode_line_analysis.sh` ]

- `combined/plots/Heatmap_zscore_<samplefile_name>_on_<regionfile_name>.pdf`\
Clustered heatmap of all the differentially expressed genes between all pairs of RNAseq samples in the `<samplefile_name>_analysis_samplefile.txt`, scaling each row by zscore among the RNAseq replicate samples (highlights differences between samples).

- `combined/plots/<samplefile_name>_<regionfile_name>_heatmap_regions.pdf`\
Heatmap of the enrichment for all ChIP and RNAseq samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, scaling each region to the same length, in decreasing order of overall enrichment in all samples.\
Tool details: https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html \
\[ generated in R by `MaizeCode_R_DEG.r`, started from `MaizeCode_line_analysis.sh` ]

- `combined/plots/<samplefile_name>_on_<regionfile_name>_heatmap_regions_k5.pdf`\
Heatmap of the enrichment for all ChIP and RNAseq samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, scaling each region to the same length, clustered into 5 regions by kmeans.\
Tool details: https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html \
\[ generated in R by deeptools from `MaizeCode_line_analysis.sh` ]

- `combined/plots/<samplefile_name>_on_<regionfile_name>_heatmap_tss.pdf`\
Heatmap of the enrichment for all ChIP and RNAseq samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, aligning all regions by their transcription start site, in decreasing order of overall enrichment in all samples.\
Tool details: https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html \
\[ generated in R by deeptools from `MaizeCode_line_analysis.sh` ]

- `combined/plots/<samplefile_name>_on_<regionfile_name>_heatmap_tss_k5.pdf`\
Heatmap of the enrichment for all ChIP and RNAseq samples in the `<samplefile_name>_analysis_samplefile.txt` on the regions from `<regionfile_name>.bed`, aligning all regions by their transcription start site, clustered into 5 regions by kmeans.\
Tool details: https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html \
\[ generated in R by deeptools from `MaizeCode_line_analysis.sh` ]

- `combined/plots/<samplefile_name>_on_<regionfile_name>_heatmap_DEG.pdf`\
Heatmap of the enrichment for all ChIP samples in the `<samplefile_name>_analysis_samplefile.txt` on the groups of differentially expressed genes called between the all pairs of RNAseq samples, scaling each gene to the same length, in decreasing order of overall enrichment in each group of UP and DOWN-regulated genes, using a specific scale for each mark (warning: genes differentially expressed between several pairs of tissues will be present in the corresponding clusters).\
Tool details: https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html \
\[ generated in R by deeptools from `MaizeCode_line_analysis.sh` ]

- `combined/plots/<samplefile_name>_on_<regionfile_name>_heatmap_all_DEGs.pdf`\
Heatmap of the enrichment for all ChIP samples in the `<samplefile_name>_analysis_samplefile.txt` on all the differentially expressed genes called between the all pair of RNAseq samples, scaling each gene to the same length, and clustering into 5 groups by kmeans, using a specific scale for each mark (warning: the generated clusters of genes are not linked to the samples they were originally called in as differentially expressed).\
Tool details: https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html \
\[ generated in R by deeptools from `MaizeCode_line_analysis.sh` ]


**-l tmp_free to be readded after Todd works on it**
line analysis = 100G
maizecode = 10G
maizecode_analysis = 1G


