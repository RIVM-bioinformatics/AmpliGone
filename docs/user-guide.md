# User Guide

AmpliGone works with both reads in FASTQ format, as well as aligned data in BAM-format. However, when data is presented in the BAM-format then only read-data (sequence and quality scores) will be used. Other data present in the BAM-format will not be used in this version of AmpliGone.

Currently, AmpliGone supports Nanopore data and Illumina data. The Illumina platform (NextSeq/MiSeq/HiSeq/other) does not matter.

It is however important that you know the read-length in relation to the amplicon length.
AmpliGone expects this information in the form of an 'amplicon-type'.
Please see [this page](amplicon-types.md) for more information regarding supported amplicon-types.

## Preparing input data

For optimal results we ask you to format your inputs to make sure you get the best results in your analysis.

It is required that adapters/barcodes have been removed from your sequencing reads *before* they are processed by AmpliGone.
Processing read-data with adapters/barcodes still attached to your reads will result in inaccurate output data.

Additionally, AmpliGone works best with reads that have already been processed by quality filtering tools such as [fastp](https://github.com/OpenGene/fastp) or [Trimmomatic](https://github.com/usadellab/Trimmomatic).


### Primers

AmpliGone searches for the primer coordinates based on the given reference. It does so by matching a given primer to the reference, up to 3 differences (substitutions) between primer and reference are tolerated before the coordinates can no longer be found.
It's therefore important that the given primers adequately match the given reference, otherwise the primer-coordinates cannot be determined.

If necessary, you can adjust the maximum amount of differences (substitutions) which are tolerated to be either more lenient or stringent. You can do so by giving the `--error-rate`/ `-er` flag in your command, followed by a single number. For example: `--error-rate 1`

AmpliGone determines whether a given sequence is considered to be a 'forward primer' or a 'reverse primer'. This information is taken from the name of a primer sequence in the FASTA header with specific keywords, the usable keywords are as follows:

**Forward primers**

* "LEFT"
* "PLUS"
* "POSITIVE"
* "FORWARD"

**Reverse primers**

* "RIGHT"
* "MINUS"
* "NEGATIVE"
* "REVERSE"

!!! example
    ```
    >primer_1_LEFT
    ACTGGC
    >primer_2_RIGHT
    GATTCA
    ```

The primer sequences may contain IUPAC ambiguity nucleotides.

### Reference

Unlike the primer sequences, AmpliGone currently only supports a reference that *do not* contain IUPAC ambiguity nucleotides. This may change in a future version.

## Basic usage example

**Example 1**: removing primers from Nanopore "end-to-end"-amplicon reads, using 8 threads:

```bash
ampligone \
    --input input.fastq \
    --output output.fastq \
    --reference reference.fasta \
    --primers primers.fasta \
    --amplicon-type end-to-end \
    --threads 8
```

**Example 2**: removing primers from Illumina "end-to-mid"-amplicon reads, using 12 threads:

```bash
ampligone \
    --input input.fastq \
    --output output.fastq \
    --reference reference.fasta \
    --primers primers.fasta \
    --amplicon-type end-to-mid \
    --threads 12
```

## Exporting found primer coordinates

In a downstream analysis you might want to know which primers have actually been removed from which coordinates by AmpliGone.
AmpliGone can provide you with this information with the `--export-primers {file}` flag, replace `{file}` with your desired output file.
Using this flag will give you a BED-file with found primer coordinates as in the example below:

| ref | start | end | name | score | strand |
| ---- | ---- | ---- | ---- | ---- | ---- |
NC_045512.2 | 30 | 54 | primer_1_LEFT | . | + |
NC_045512.2 | 385 | 410 | primer_1_RIGHT | . | - |
NC_045512.2 | 320 | 342 | primer_2_LEFT | . | + |
NC_045512.2 | 704 | 726 | primer_2_RIGHT | . | - |

!!! example "Example command"
    ```bash
    ampligone \
        --input input.fastq \
        --output output.fastq \
        --reference reference.fasta \
        --primers primers.fasta \
        --amplicon-type end-to-end \
        --export-primers removed_coordinates.bed \
        --threads 12
    ```

## Using multiple threads

AmpliGone is multi-threaded by design, as it splits your input data to process over the amount of threads which are available for use.
Therefore, AmpliGone scales (almost) linearly when given more threads. A computer with a lot of computing power is therefore advised.

AmpliGone defaults to the amount of threads which are available in your system. (if the CPU of your computer has 24 threads, AmpliGone will use all 24 threads by default)
You can use the `--threads` or `-t` flag to set a different number of threads to use.


## Using AmpliGone in a pipeline/workflow

In pipelines it's often required to that output-files are always created, even if they are empty.
AmpliGone has the `-to` flag for this. Which stands for 'touch-outputs'.

When this flag is given, AmpliGone will *always* create the output-files. Even if AmpliGone would normally fail because of (for example) an empty input file or when primer-sequences couldn't be found on the reference.

!!! warning "Please use with caution"
    Using the `-to` flag will ensure that output-files are always created.
    But because of this, common errors may go unnoticed more easily. Please test your experiment setup manually before integrating it in a pipeline
