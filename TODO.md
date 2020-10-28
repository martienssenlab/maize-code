# To do list:

- [ ] Update README file
- [x] Add details on the input files (gff3 and fasta)
- [x] Rename deeptools folder
- [x] Change type variable to datatype
- [ ] Change hard-linked paths to relative for scripts (alternative to dirname "$0")
- [ ] Work on intermediate file hygiene
- [ ] Work on the compatiblity with other systems (no cluster, no qsub ...)
- [ ] Compute PBC1 and PBC2 (PCR Bottlenecking Coefficient) (https://www.encodeproject.org/data-standards/terms/#library)
- [ ] Make genome coverage files plots for ChIPseq samples
- [ ] Improve portability/management of R packages
- [ ] Manage the way to get the fastq files depending on where they are coming from (server, sra, ...)
- [x] Add the references in the sample file to be able to map to different references form the same samplefile and can then be processed together later on +
get the data type from the samplefile too (if H* or Input it is ChIP, if *RNA* or RAMPAGE it is RNA)
- [ ] Make the RNA pipeline
- [ ] Continue analysis pipeline
- [ ] Decide on changing naming of files to record reference they have been mapped to? (for the moment, if a file has been mapped to reference A, it will not be mapped on A again, but if asked to be mapped on reference B, new files will overwrite existing reference A-mapped files!)
- [ ] Fix warning when asked to check for the existence of a file that has several possiblities (e.g. [ ! -e ./$datatype/fastq/${name}*.fastq.gz ] for PE data in MaizeCode.sh)
- [ ] Add output information in the readme (files, logs, summaries, plots, etc..)
