To do list:

0) Update README file

1) DONE Add details on the input files (gff3 and fasta)

2) DONE (Rename deeptools folder)

3) DONE (Change type variable to datatype)

4) Change hard-linked paths to relative, most importantly for scripts (Evan)

5) Intermediate file hygiene

6) Down the line? Compatiblity for potential user that has a different set-up (no cluster, no qsub ...)

7) DONE Add the references in the sample file to be able to map to different references form the same samplefile and can then be processed together later on +
get the data type from the samplefile too (if H* or Input it is ChIP, if *RNA* or RAMPAGE it is RNA)

8) Compute PBC1 and PBC2 (PCR Bottlenecking Coefficient) (https://www.encodeproject.org/data-standards/terms/#library)

9) Management of R packages

10) Manage the way to get the fastq files depending on where they are coming from (server, sra, ...)

11) Make the RNA mapping+QC pipeline

12) Make a wrapper script to incorporate RNA samples into the analysis
