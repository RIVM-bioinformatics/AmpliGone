import sys

from setuptools import find_packages, setup

from AmpliGone.version import __version__

if sys.version_info.major != 3 or sys.version_info.minor < 7:
    print("Error: you must execute setup.py using Python 3.7 or later")
    sys.exit(1)

with open("README.md", "rb") as readme:
    DESCR = readme.read().decode()


setup(
    name="AmpliGone",
    version=__version__,
    author="Florian Zwagemaker",
    author_email="ids-bioinformatics@rivm.nl",
    license="AGPLv3",
    packages=find_packages(),
    install_requires=[
        "pysam>=0.16",
        "pandas>=1.2.3",
        "numpy>=1.20",
        "mappy>=2.17",
        "biopython>=1.78",
        "parmap>=1.5.2",
        "tqdm>=4.59.0",
    ],
    entry_points={
        "console_scripts": [
            "ampligone = AmpliGone.AmpliGone:main",
            "AmpliGone = AmpliGone.AmpliGone:main",
        ]
    },
    keywords=[],
    zip_safe=False,
)
