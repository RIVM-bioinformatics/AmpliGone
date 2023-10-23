# Installation instructions

AmpliGone is only available on Linux (or Linux-based) operating systems. MacOS may also work but is not tested. AmpliGone will *not work* on Windows.

AmpliGone can be installed easily with [Conda](https://anaconda.org/bioconda/ampligone) or [pip](https://pypi.org/project/AmpliGone/).

### Installation with conda
AmpliGone releases are published in Bioconda, and can be installed with the following command:


```bash
conda install -c bioconda ampligone
```

!!! tip
    You can also use `mamba` as a replacement for installing conda packages, use the following command if you have `mamba` installed on your system
    ```bash
    mamba install -c bioconda ampligone
    ```

Installation through conda/mamba is the recommended method for installing AmpliGone

### Installation with pip
AmpliGone releases are also published on [PyPI](https://pypi.org/project/AmpliGone/), and can be installed using pip using the following command:

```bash
pip install AmpliGone
```

### Download and installation from source

!!! tip "Make a new Conda environment before continuing"
    If you have Conda installed on your system, please create and activate a new environment before continuing.

    Use the following command to create and activate a new Conda environment named "AmpliGone" based on the environment-recipe we provide in the github-repository

    ```bash
    conda env create -f env.yml; conda activate AmpliGone
    ```

    The "AmpliGone" conda-environment should now be active

First start by cloning the repository and make sure you're on the latest released version of AmpliGone, then install the downloaded package with pip:

```bash
# clone the repository from github and switch to the latest released version
git clone https://github.com/RIVM-bioinformatics/AmpliGone.git; cd AmpliGone; git checkout tags/$(git tag --sort=committerdate | tail -1) >> /dev/null
# Install the downloaded package with pip
pip install .
```

AmpliGone should now be installed!  
You can verify if installation was successful by typing `ampligone --version` on the command-line, this should show the installed AmpliGone version.

## Pipeline/workflow integration

You can easily integrate AmpliGone in your SnakeMake or Nextflow bioinformatics workflow if you use Conda environments in your workflow.

To do this, add `ampligone` and the correct channels to your conda-environment recipe as below:

```yaml
name: example-env
channels:
    - bioconda
    - conda-forge
dependencies:
    - ampligone
```

Conda will now install AmpliGone and its dependencies in the specified snakemake conda-environment.
