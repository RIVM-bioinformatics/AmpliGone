import sys

from setuptools import find_packages, setup

from AmpliGone import __prog__, __version__

if sys.version_info.major != 3 or sys.version_info.minor < 10:
    print("Error: you must execute setup.py using Python 3.10 or later")
    sys.exit(1)

with open("README.md", "r", encoding="utf-8") as readme:
    DESCR = readme.read()

setup(
    name=__prog__,
    version=__version__,
    url="https://rivm-bioinformatics.github.io/AmpliGone/",
    project_urls={"Source Code": "https://github.com/RIVM-bioinformatics/AmpliGone"},
    author="Florian Zwagemaker",
    author_email="ids-bioinformatics@rivm.nl",
    description=f"{__prog__} is a tool which accurately removes primer sequences from FastQ NGS reads in an amplicon sequencing experiment",
    long_description=DESCR,
    long_description_content_type="text/markdown",
    python_requires=">=3.10",
    license="AGPLv3",
    packages=find_packages(),
    install_requires=[
        "pysam==0.22.*",
        "pandas==2.2.*",
        "mappy==2.28",
        "biopython==1.84",
        "parmap==1.7.*",
        "parasail==1.3.4",
        "rich==13.7.*",
        "pgzip==0.3.4",
    ],
    entry_points={
        "console_scripts": [
            "ampligone = AmpliGone.__main__:main",
            "AmpliGone = AmpliGone.__main__:main",
        ]
    },
    keywords=[],
    zip_safe=False,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
