import sys

from setuptools import find_packages, setup

from AmpliGone.version import __version__

if sys.version_info.major != 3 or sys.version_info.minor < 8:
    print("Error: you must execute setup.py using Python 3.8 or later")
    sys.exit(1)

with open("README.md", "r", encoding="utf-8") as readme:
    DESCR = readme.read()

setup(
    name="AmpliGone",
    version=__version__,
    url="https://rivm-bioinformatics.github.io/AmpliGone/",
    project_urls={"Source Code": "https://github.com/RIVM-bioinformatics/AmpliGone"},
    author="Florian Zwagemaker",
    author_email="ids-bioinformatics@rivm.nl",
    description="Ampligone is a tool which accurately removes primer sequences from FastQ NGS reads in an amplicon sequencing experiment",
    long_description=DESCR,
    long_description_content_type="text/markdown",
    python_requires=">=3.8",
    license="AGPLv3",
    packages=find_packages(),
    install_requires=[
        "pysam==0.19.*",
        "pandas>=1.3.0",
        "mappy==2.24",
        "biopython>=1.79",
        "parmap==1.5.*",
        "regex>=2021.11.10",
        "rich==12.5.*",
    ],
    entry_points={
        "console_scripts": [
            "ampligone = AmpliGone.AmpliGone:main",
            "AmpliGone = AmpliGone.AmpliGone:main",
        ]
    },
    keywords=[],
    zip_safe=False,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
