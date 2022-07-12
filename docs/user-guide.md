# User Guide

AmpliGone works with both reads in FASTQ format, as well as aligned data in BAM-format. However, when data is presented in the BAM-format then only read-data (sequence and quality scores) will be used. Other data present in the BAM-format will not be used in this version of AmpliGone.

Currently, AmpliGone supports Nanopore data and Illumina data. The Illumina platform (NextSeq/MiSeq/HiSeq/other) does not matter.

It is however important that you know the read-length in relation to the amplicon length.
AmpliGone expects this information in the form of an 'amplicon-type'.
Please see [this page](amplicon-types.md) for more information regarding supported amplicon-types.

## Preparing input data

For optimal results we ask you to format your inputs to make sure you get the best results in your analysis.

It is required that adapters/barcodes have been removed from your sequencing reads *before* they are processed by AmpliGone.
Processing read-data with adapters/barcodes still attached to your reads may result in inaccurate output data.

Additionally, AmpliGone works best with reads that have already been processed by quality filtering tools such as [fastp](https://github.com/OpenGene/fastp) or [Trimmomatic](https://github.com/usadellab/Trimmomatic).


### Primers

AmpliGone has two options for primer input: [BED](https://en.wikipedia.org/wiki/BED_(file_format)) and [fasta](https://en.wikipedia.org/wiki/FASTA_format).
The BED format specifies the coordinates of the primers with respect to the given reference genome. If this is the format you use, you can skip to [Primer Orientation](#primer-orientation)

When the primers are supplied in fasta format, AmpliGone searches for the primer coordinates based on the given reference.
By default, mismatches are tolerated for up to 10% (0.1) of the length of the primer sequence.
It's therefore important that the given primers adequately match the given reference, otherwise the primer-coordinates cannot be determined.

If necessary, you can adjust the maximum amount of differences (substitutions) which are tolerated to be either more lenient or stringent. You can do so by giving the `--error-rate`/ `-er` flag in your command, followed by a decimal number. For example `--error-rate 0.15` can be used to make the search more lenient.  
Please see the [usage examples](#basic-usage-examples) to see this flag in combination with other settings for AmpliGone.

### Primer Orientation
AmpliGone determines whether a given sequence is considered to be a 'forward primer' or a 'reverse primer'. This information is taken from the name of a primer sequence in the FASTA header with specific keywords, the usable keywords are as follows:

=== "Usable keywords"
    | Forward primers | Reverse primers |
    | --------------- | --------------- |
    | LEFT            | RIGHT           | 
    | PLUS            | MINUS           |
    | POSITIVE        | NEGATIVE        |
    | FORWARD         | REVERSE         | 

=== "Example fasta" 
    ```
    >primer_1_LEFT
    ACTGGC
    >primer_1_RIGHT
    TGGCTCA
    >primer_2_FORWARD
    ACAATTCG
    >primer_2_REVERSE
    TATTAAGC
    ```

The primer sequences may contain IUPAC ambiguity nucleotides, though be aware that this may result in differently found coordinates as expected depending on your experiments.

### Reference

Unlike the primer sequences, AmpliGone currently only supports a reference sequence that *<u>does not</u>* contain IUPAC ambiguity nucleotides. This may change in a future version.

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

Please see the [usage examples](#basic-usage-examples) to see this flag in combination with other settings for AmpliGone.

## Using multiple threads

AmpliGone is multi-threaded by design, as it splits your input data to process over the amount of threads which are available for use.
Therefore, AmpliGone scales (almost) linearly when given more threads. A computer with a lot of computing power is therefore advised.

AmpliGone defaults to the amount of threads which are available in your system. (if the CPU of your computer has 24 threads, AmpliGone will use all 24 threads by default)
You can use the `--threads` or `-t` flag to set a different number of threads to use.

## Basic usage examples

=== "End-to-end amplicons"

    Here you can find various examples to remove primers from **"end-to-end"** amplicon data.  
    Please check what an "end-to-end" amplicon is <u>[here](amplicon-types.md#end-to-end)</u> before continuing.

    ???+ Example "Basic primer removal with a Fasta-file as primer-input"

        * The primers in this command are provided in **Fasta** format.  
        * The default setting for a primer-mismatches (10%) will be used as no other value is provided in this command.  
        * Aside from basic logging messages, AmpliGone will not output any additional data such as removed primer coordinates.  
        * The flag `--threads` is *not* provided in this command, AmpliGone will therefore use *all* available threads on the system.
        
        ```bash
        ampligone \
            --input input.fastq \
            --output output.fastq \
            --reference reference.fasta \
            --primers primers.fasta \
            --amplicon-type end-to-end
        ```

    ??? Example "Basic primer removal with a BED-file as primer-input"

        * The primers in this command are provided in **BED** format.  
        * The default setting for a primer-mismatches (10%) will be used as no other value is provided in this command.  
        * Aside from basic logging messages, AmpliGone will not output any additional data such as removed primer coordinates.
        * The flag `--threads 6` is given in this command, AmpliGone will therefore use only 6 threads.

        ```bash
        ampligone \
            --input input.fastq \
            --output output.fastq \
            --reference reference.fasta \
            --primers primer_information.bed \
            --amplicon-type end-to-end \
            --threads 6
        ```
    
    ??? Example "Primer removal with a BED-file as primer-input and exporting removed primers"

        * The primers in this command are provided in **BED** format.
        * The default setting for a primer-mismatches (10%) will be used as no other value is provided in this command.  
        * AmpliGone will output both basic logging messages as an extra output file containing the found and removed primer coordinates in **BED** format.
        * The flag `--threads 8` is given in this command, AmpliGone will therefore use only 8 threads.

        ```bash
        ampligone \
            --input input.fastq \
            --output output.fastq \
            --reference reference.fasta \
            --primers primers.fasta \
            --amplicon-type end-to-end \
            --export-primers removed_primers.bed \
            --threads 8
        ```
    
    ??? Example "Primer removal with a Fasta-file as primer-input, using a custom primer mismatch rate and exporting removed primers"

        * The primers in this command are provided in **Fasta** format.
        * The flag `--error-rate 0.2` is given in this command, meaning that a primer mismatch rate of 20% will be tolerated by AmpliGone instead of the default 10%
        * AmpliGone will output both basic logging messages as an extra output file containing the found and removed primer coordinates in **BED** format.
        * The flag `--threads 16` is given in this command, AmpliGone will therefore use only 16 threads.

        ```bash
        ampligone \
            --input input.fastq \
            --output output.fastq \
            --reference reference.fasta \
            --primers primers.fasta \
            --amplicon-type end-to-end \
            --error-rate 0.2 \
            --export-primers removed_primers.bed \
            --threads 16
        ```

=== "End-to-mid amplicons"

    Here you can find various examples to remove primers from **"end-to-mid"** amplicon data.  
    Please check what an "end-to-mid" amplicon is <u>[here](amplicon-types.md#end-to-mid)</u> before continuing.

    ???+ Example "Basic primer removal with a Fasta-file as primer-input"

        * The primers in this command are provided in **Fasta** format.  
        * The default setting for a primer-mismatches (10%) will be used as no other value is provided in this command.  
        * Aside from basic logging messages, AmpliGone will not output any additional data such as removed primer coordinates.  
        * The flag `--threads` is *not* provided in this command, AmpliGone will therefore use *all* available threads on the system.

        ```bash
        ampligone \
            --input input.fastq \
            --output output.fastq \
            --reference reference.fasta \
            --primers primers.fasta \
            --amplicon-type end-to-mid
        ```

    ??? Example "Basic primer removal with a BED-file as primer-input"

        * The primers in this command are provided in **BED** format.  
        * The default setting for a primer-mismatches (10%) will be used as no other value is provided in this command.  
        * Aside from basic logging messages, AmpliGone will not output any additional data such as removed primer coordinates.
        * The flag `--threads 6` is given in this command, AmpliGone will therefore use only 6 threads.

        ```bash
        ampligone \
            --input input.fastq \
            --output output.fastq \
            --reference reference.fasta \
            --primers primer_information.bed \
            --amplicon-type end-to-mid \
            --threads 6
        ```

    ??? Example "Primer removal with a BED-file as primer-input and exporting removed primers"

        * The primers in this command are provided in **BED** format.
        * The default setting for a primer-mismatches (10%) will be used as no other value is provided in this command.  
        * AmpliGone will output both basic logging messages as an extra output file containing the found and removed primer coordinates in **BED** format.
        * The flag `--threads 8` is given in this command, AmpliGone will therefore use only 8 threads.

        ```bash
        ampligone \
            --input input.fastq \
            --output output.fastq \
            --reference reference.fasta \
            --primers primers.fasta \
            --amplicon-type end-to-mid \
            --export-primers removed_primers.bed \
            --threads 8
        ```
    
    ??? Example "Primer removal with a Fasta-file as primer-input, using a custom primer mismatch rate and exporting removed primers"

        * The primers in this command are provided in **Fasta** format.
        * The flag `--error-rate 0.2` is given in this command, meaning that a primer mismatch rate of 20% will be tolerated by AmpliGone instead of the default 10%
        * AmpliGone will output both basic logging messages as an extra output file containing the found and removed primer coordinates in **BED** format.
        * The flag `--threads 16` is given in this command, AmpliGone will therefore use only 16 threads.

        ```bash
        ampligone \
            --input input.fastq \
            --output output.fastq \
            --reference reference.fasta \
            --primers primers.fasta \
            --amplicon-type end-to-mid \
            --error-rate 0.2 \
            --export-primers removed_primers.bed \
            --threads 16
        ```

## Using AmpliGone in a pipeline/workflow

In pipelines it's often required to that output-files are always created, even if they are empty.
AmpliGone has the `-to` flag for this. Which stands for 'touch-outputs'.

When this flag is given, AmpliGone will *always* create the output-files. Even if AmpliGone would normally fail because of (for example) an empty input file or when primer-sequences couldn't be found on the reference.

!!! warning "Please use with caution"
    Using the `-to` flag will ensure that output-files are always created.
    But because of this, common errors may go unnoticed more easily. Please test your experiment setup manually before integrating it in a pipeline
