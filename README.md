# AmpliGone

AmpliGone finds and removes primer sequences from PCR amplicons based on coordinates.  
AmpliGone works for Illumina and Nanopore data and supports both 'end-to-end' as well as 'end-to-mid' amplicons to be cleaned up.

## Usage:

* Example 1 - Using AmpliGone to remove primers from Nanopore end-to-end amplicons:
```
ampligone -i <FASTQ FILE> -o <OUTPUT FASTQ FILE> -ref <REFERENCE FASTA> -pr <PRIMER FASTA> --amplicon-type end-to-end --threads 12
```

* Example 2 - Using AmpliGone to remove primers from Illumina end-to-mid amplicons:
```
ampligone -i <FASTQ FILE> -o <OUTPUT FASTQ FILE> -ref <REFERENCE FASTA> -pr <PRIMER FASTA> --amplicon-type end-to-mid --threads 12
```