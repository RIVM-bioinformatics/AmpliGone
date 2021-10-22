---
hide:
  - navigation
  - toc
---

[![CodeFactor](https://www.codefactor.io/repository/github/rivm-bioinformatics/ampligone/badge)](https://www.codefactor.io/repository/github/rivm-bioinformatics/ampligone)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/RIVM-bioinformatics/AmpliGone?include_prereleases)
![GitHub](https://img.shields.io/github/license/RIVM-bioinformatics/AmpliGone)
# AmpliGone

AmpliGone is a tool which accurately finds and removes primer sequences from NGS reads in an amplicon experiment.

In contrast to a lot of other primer-removal tools, AmpliGone does not actively look for primer sequences within the NGS reads. Instead, reads are trimmed based on primer sequence coordinates in relation to a given reference sequence.  
Additionally, AmpliGone is able to compensate for, and therefore properly clean, reads that start or end outside of a primer-region as this is a common occurrence in amplicon-based sequencing data.

AmpliGone is build and tested with Nanopore and Illumina data (fastq) in mind and supports both 'end-to-end' as well as 'end-to-mid' amplicons to be cleaned. Please see [this page](amplicon-types.md) to learn more about this terminology.

AmpliGone is available under the [AGPLv3 licence](https://www.gnu.org/licenses/agpl-3.0.en.html) 