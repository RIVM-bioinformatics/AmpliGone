# AmpliGone

AmpliGone is a tool which accurately finds and removes primer sequences from NGS reads in an amplicon experiment.

In contrast to a lot of other primer-removal tools, AmpliGone does not actively look for primer sequences within the NGS reads. Instead, reads are trimmed based on primer sequence coordinates in relation to a given reference sequence.
Additionally, AmpliGone is able to compensate for, and therefore properly clean, reads that start or end outside of a primer-region as this is a common occurrence in amplicon-based sequencing data.

AmpliGone is build and tested with Nanopore and Illumina data (fastq) in mind and supports both 'end-to-end' as well as 'end-to-mid' amplicons to be cleaned.  
Please see [this page](https://rivm-bioinformatics.github.io/AmpliGone/latest/amplicon-types/) to learn more about this terminology.

Please see the [documentation](https://rivm-bioinformatics.github.io/AmpliGone/) for more information as well as [installation instructions](https://rivm-bioinformatics.github.io/AmpliGone/latest/installation/) and [usage instructions](https://rivm-bioinformatics.github.io/AmpliGone/latest/user-guide/).

AmpliGone is available under the [AGPLv3 license](https://www.gnu.org/licenses/agpl-3.0.en.html).