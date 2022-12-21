[![CodeFactor](https://www.codefactor.io/repository/github/rivm-bioinformatics/ampligone/badge)](https://www.codefactor.io/repository/github/rivm-bioinformatics/ampligone)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/RIVM-bioinformatics/AmpliGone?include_prereleases)
![GitHub](https://img.shields.io/github/license/RIVM-bioinformatics/AmpliGone)  
[![Build and release](https://github.com/RIVM-bioinformatics/AmpliGone/actions/workflows/release.yml/badge.svg)](https://github.com/RIVM-bioinformatics/AmpliGone/actions/workflows/release.yml)
![GitHub deployments](https://img.shields.io/github/deployments/RIVM-bioinformatics/AmpliGone/github-pages?label=Documentation%20deployment)  
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ampligone/README.html)


# AmpliGone

AmpliGone is a tool which accurately finds and removes primer sequences from NGS reads in an amplicon experiment.

In contrast to a lot of other primer-removal tools, AmpliGone does not actively look for primer sequences within the NGS reads. Instead, reads are trimmed based on primer sequence coordinates in relation to a given reference sequence.
Additionally, AmpliGone is able to compensate for, and therefore properly clean, reads that start or end outside of a primer-region as this is a common occurrence in amplicon-based sequencing data.

AmpliGone is build and tested with Nanopore and Illumina data (fastq) in mind and supports 'end-to-end', 'end-to-mid' and 'fragmented' amplicons to be cleaned.  
Please see [this page](https://rivm-bioinformatics.github.io/AmpliGone/latest/amplicon-types/) to learn more about this terminology.

## Installation instructions

AmpliGone can be installed easily with [conda](https://anaconda.org/bioconda/ampligone) or [pip](https://pypi.org/project/AmpliGone/).
Installation through conda is recommended.

### Installation with conda

```bash
conda install -c bioconda ampligone
```

### Installation with pip

```bash
pip install AmpliGone
```

Please see the [documentation](https://rivm-bioinformatics.github.io/AmpliGone/) for more information as well as extended [installation instructions](https://rivm-bioinformatics.github.io/AmpliGone/latest/installation/) and [usage instructions](https://rivm-bioinformatics.github.io/AmpliGone/latest/user-guide/).

AmpliGone is freely available under the [AGPLv3 license](https://www.gnu.org/licenses/agpl-3.0.en.html).