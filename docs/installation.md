# Installation instructions

AmpliGone is only available on Linux (or Linux-based) operating systems. MacOS may also work but is not tested. AmpliGone will *not work* on Windows.

AmpliGone is available via [PyPI](https://pypi.org/project/AmpliGone/), and and can be installed using pip:

```bash
pip install AmpliGone
```

???+ note
    AmpliGone is not yet available via Conda. We will update these docs when this is available.

## Prerequisites

AmpliGone requires Python 3.8 or later to be installed on your system (or in an environment).

Other dependencies will be installed during setup, you don't have to install them manually. These extra dependencies are as follows:

* pysam ≥ 0.16
* pandas ≥ 1.2.3
* numpy ≥ 1.20
* mappy ≥ 2.17
* biopython ≥ 1.79
* parmap ≥ 1.5.2
* regex ≥ 2021.11.10
* tqdm ≥ 4.59.0

We strongly advise you to use a conda environment (or similar) to make sure there won't be any conflicts in package dependencies.

## Download and install from source

First start by cloning the repository and make sure you're on the latest released version of AmpliGone:
```bash
git clone https://github.com/RIVM-bioinformatics/AmpliGone.git; cd AmpliGone; git checkout tags/$(git tag --sort=committerdate | tail -1) >> /dev/null
```

You're now in the newly created "AmpliGone" directory.

!!! tip "Make a new Conda environment before continuing"
    If you have Conda installed on your system, please create and activate a new environment before continuing.

    Use the following command to create and activate a new Conda environment named "AmpliGone" based on the environment-recipe we provide in the github-repository

    ```bash
    conda env create -f env.yml; conda activate AmpliGone
    ```

    The "AmpliGone" conda-environment should now be active 

You can now install AmpliGone via the following command:
```bash
pip install .
```

AmpliGone should now be installed!  
You can verify if installation was successful by typing `ampligone --version` on the command-line, this should show the installed AmpliGone version.

## Pipeline/workflow integration

You can easily integrate AmpliGone in your snakemake bioinformatics workflow if you use Conda environments in your workflow.

To do this, simply add the following structure to your conda-environment recipe, replace `{VERSION}` with the AmpliGone version you wish to use:

```yaml
dependencies:
    - pip
    - pip:
        - AmpliGone
```

Conda will now install AmpliGone and its dependencies in the specified snakemake conda-environment.
