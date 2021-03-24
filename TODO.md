# To do list:

- [ ] Update README file
- [x] Add details on the input files (gff3 and fasta)
- [x] Rename deeptools folder
- [x] Change type variable to datatype
- [x] Change hard-linked paths to relative for scripts (alternative to dirname "$0")
- [x] Store reference genomes on server for normalization and make a script to download them if needed, that can be edited for downloading from the web: reference folder exists, but needs to be copied first, not sure we want to create a new genome folder everytime. Maybe can be added as an option (if genome folder exists, use it, if not, download?)
- [ ] Work on intermediate file hygiene (many intermediate files that could be removed). Maybe add option -output (or something) to chose what should be kept?
- [ ] Work on the compatiblity with other systems (no cluster, no qsub ...)
- [ ] Compute PBC1 and PBC2 (PCR Bottlenecking Coefficient) (https://www.encodeproject.org/data-standards/terms/#library)
- [ ] Make genome coverage files plots for ChIPseq samples
- [ ] Improve portability/management of R packages
- [x] Manage the way to get the fastq files depending on where they are coming from (server, sra, ...)
- [ ] Archive raw fastq files
- [x] Add the references in the sample file to be able to map to different references form the same samplefile and can then be processed together later on +
get the data type from the samplefile too (if H* or Input it is ChIP, if \*RNA* or RAMPAGE it is RNA)
- [x] Make the RNA pipeline
- [ ] Make the shRNA pipeline
- [ ] Continue analysis pipeline
- [ ] Decide on changing naming of files to record reference they have been mapped to? (for the moment, if a file has been mapped to reference A, it will not be mapped on A again, but if asked to be mapped on reference B, new files will overwrite existing reference A-mapped files!)
- [x] Fix warning when asked to check for the existence of a file that has several possiblities (e.g. [ ! -e ./$datatype/fastq/${name}*.fastq.gz ] for PE data in MaizeCode.sh)
- [x] Add output information in the readme (files, logs, summaries, plots, etc..)
- [ ] Improve logs (callback if errors, follow-up on where it's actually at, naming, etc...). The 'wait' command will break if there is an error instead of continue with checking the touch files.
- [x] Decide how to work with ChIP replicates For many downstream analysis, it would be much easier to have one file for each line * tissue * mark. Do we do it by default on peaks called after merging bam files, on peaks passing the IDR threshold (very limited when samples are not great), or peaks that are called in the _best_ replicate (to be defined). This needs to be automatize because going through all the different combinations possible is overwhelming. __Chose pseudo-replicates method__
- [ ] Make a genome browser session? (trackhub for UCSC? Cyverse?)
- [ ] Improve flexibility of scripts based on number of replicates
- [ ] Increase efficiency of the whole pipeline? As of now, the whole pipeline on all ChIP samples for B73 roots took ~12h to run.
- [ ] Add "CC" ontology to GO analysis (weird NA error) and work around corner cases for scatter plots
