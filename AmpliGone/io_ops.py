import gzip
import pathlib
import sys
from typing import Any, Dict, Hashable, List, TextIO, Tuple

import pandas as pd
import pysam

from .func import log


def is_zipped(filename: str) -> bool:
    """
    Check if the given file is a gzipped file.

    Parameters
    ----------
    filename : str
        The name of the file to be checked.

    Returns
    -------
    bool
        True if the file is gzipped, False otherwise.

    Examples
    --------
    >>> is_zipped("file.txt.gz")
    True

    >>> is_zipped("file.txt")
    False
    """
    return bool(".gz" in pathlib.Path(filename).suffixes)


def is_fastq(filename: str) -> bool:
    """
    Check if the given file is a FASTQ file.

    Parameters
    ----------
    filename : str
        The name of the file to be checked.

    Returns
    -------
    bool
        True if the file is a FASTQ file, False otherwise.

    Examples
    --------
    >>> is_fastq("file.fastq")
    True

    >>> is_fastq("file.txt")
    False
    """
    ext = [".fastq", ".fq"]
    return bool(any(item in ext for item in pathlib.Path(filename).suffixes))


def is_bam(filename: str) -> bool:
    """
    Check if the given file is a BAM file.

    Parameters
    ----------
    filename : str
        The name of the file to be checked.

    Returns
    -------
    bool
        True if the file is a BAM file, False otherwise.

    Examples
    --------
    >>> is_bam("file.bam")
    True

    >>> is_bam("file.txt")
    False
    """
    return bool(".bam" in pathlib.Path(filename).suffixes)


def read_gzip(filename: str) -> TextIO:
    """
    Open a gzip file for reading and return an opened file object.

    Parameters
    ----------
    filename : str
        The name of the file to be opened.

    Returns
    -------
    TextIO
        An opened file object.

    Examples
    --------
    >>> with read_gzip("file.txt.gz") as f:
    ...     print(f.read())
    ...
    This is a gzipped file.

    """
    return gzip.open(filename, "rt")


def read_fastq(filename: str) -> TextIO:
    """
    Open a FASTQ file for reading and return an opened file object.

    Parameters
    ----------
    filename : str
        The name of the file to be opened.

    Returns
    -------
    TextIO
        An opened file object.

    Examples
    --------
    >>> with read_fastq("file.fastq") as f:
    ...     print(f.read())
    ...
    @SEQ_ID
    GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    +
    !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

    """
    return open(filename, "rt")


def fastq_opener(inputfile: str) -> TextIO:
    """
    Open a FASTQ file for reading and return an opened file object.
    If the file is gzipped, it will be unzipped before opening.

    Parameters
    ----------
    inputfile : str
        The name of the input file.

    Returns
    -------
    TextIO
        An opened file object.

    Examples
    --------
    >>> with fastq_opener("file.fastq.gz") as f:
    ...     print(f.read())
    ...
    @SEQ_ID
    GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    +
    !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

    """
    if is_zipped(inputfile) is True:
        return read_gzip(inputfile)
    return read_fastq(inputfile)


def LoadBam(inputfile: str) -> pysam.AlignmentFile:
    """
    Load a BAM file and return a pysam.AlignmentFile object.

    Parameters
    ----------
    inputfile : str
        The path to the BAM file.

    Returns
    -------
    pysam.AlignmentFile
        A file object for reading the BAM file.

    Examples
    --------
    >>> bam_file = LoadBam('path/to/file.bam')
    >>> for read in bam_file:
    ...     print(read)

    """
    return pysam.AlignmentFile(inputfile, "rb")


def read_bed(filename: str) -> pd.DataFrame:
    """
    Read a BED file and return a pandas DataFrame.

    Parameters
    ----------
    filename : str
        The path to the BED file.

    Returns
    -------
    pd.DataFrame
        A DataFrame with the columns 'ref', 'start', 'end', 'name', 'score', and 'strand'.

    See Also
    --------
    pandas.read_csv : Function for reading CSV files into a DataFrame.

    Notes
    -----
    This function removes the original BED file header lines and filters out any rows that start with 'browser ' or 'track '.

    Examples
    --------
    >>> bed_df = read_bed('path/to/file.bed')
    >>> print(bed_df.head())

    """
    primer_df = pd.read_csv(
        filename,
        sep="\t",
        comment="#",
        usecols=range(6),
        header=None,
        names=["ref", "start", "end", "name", "score", "strand"],
        dtype=dict(
            ref=str,
            start="Int64",
            end="Int64",
            name=str,
            score=str,
            strand=str,
        ),
        keep_default_na=False,
    )
    primer_df = primer_df[
        ~(
            primer_df.ref.str.startswith("browser ")
            | primer_df.ref.str.startswith("track ")
        )
    ]

    return primer_df


def FlipStrand(seq: str, qual: str) -> Tuple[str, str]:
    """
    Return the reverse complement of a DNA sequence and its quality score.

    Parameters
    ----------
    seq : str
        The DNA sequence to be reverse complemented.
    qual : str
        The quality score of the DNA sequence.

    Returns
    -------
    Tuple[str, str]
        A tuple containing the reverse complement of the DNA sequence and its quality score.

    Notes
    -----
    This function uses the standard Watson-Crick base pairing rules to obtain the reverse complement of the DNA sequence.

    Examples
    --------
    >>> seq = 'ATCG'
    >>> qual = 'IIII'
    >>> flipped_seq, flipped_qual = FlipStrand(seq, qual)
    >>> print(flipped_seq, flipped_qual)
    CGAT IIII

    """
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    seq = "".join(bases)
    seq = seq[::-1]
    qual = qual[::-1]

    return seq, qual


def LoadData(inputfile: str) -> List[Tuple[str, str, str]]:
    """
    Load data from a file and return a list of tuples, where each tuple contains the name, sequence, and quality score of a read.

    Parameters
    ----------
    inputfile : str
        The input file to be processed.

    Returns
    -------
    List[Tuple[str, str, str]]
        A list of tuples. Each tuple contains the name, sequence, and quality score of a read.

    Raises
    ------
    SystemExit
        If the input file is an unsupported filetype.

    Notes
    -----
    This function supports two file types: FASTQ and BAM. If the input file is a FASTQ file, it reads the file and extracts the name, sequence, and quality score of each read. If the input file is a BAM file, it uses the LoadBam function to extract the necessary information.

    Examples
    --------
    >>> LoadData('example.fastq')
    [('read1', 'ACGT', 'IIII'), ('read2', 'TGCA', 'JJJJ')]

    >>> LoadData('example.bam')
    [('read1', 'ACGT', 'IIII'), ('read2', 'TGCA', 'JJJJ')]
    """
    Reads = []

    if is_fastq(inputfile) is True:
        with fastq_opener(inputfile) as fq:
            for line in fq:
                name: str = line.split()[0][1:]
                seq: str = next(fq).strip()
                next(fq)
                qual: str = next(fq).strip()

                Reads.append((name, seq, qual))
        return Reads
    if is_bam(inputfile) is True:
        for read in LoadBam(inputfile):
            if read.is_unmapped is True:
                continue

            name: str = read.query_name
            seq: str = read.query_sequence
            qual: str = "".join(map(lambda x: chr(x + 33), read.query_qualities))

            if read.is_reverse is True:
                seq, qual = FlipStrand(seq, qual)

            if len(seq) != len(qual):
                continue

            if len(seq) == 0:
                continue

            Reads.append((name, seq, qual))
        return Reads
    log.error(
        f'"{inputfile}" is an unsupported filetype. Please try again with a supported filetype'
    )
    sys.exit(-1)


def IndexReads(inputfile: str) -> pd.DataFrame:
    """
    Read in an input file and return a pandas dataframe with the readname, sequence, and qualities.

    Parameters
    ----------
    inputfile : str
        The path to the input file.

    Returns
    -------
    pd.DataFrame
        A pandas dataframe with the columns Readname, Sequence, and Qualities.

    Notes
    -----
    This function uses the LoadData function to extract the necessary information from the input file and returns a pandas dataframe with the extracted information.

    Examples
    --------
    >>> IndexReads('example.fastq')
           Readname Sequence Qualities
    0      read1    ACGT     IIII
    1      read2    TGCA     JJJJ
    """
    return pd.DataFrame.from_records(
        LoadData(inputfile), columns=["Readname", "Sequence", "Qualities"]
    )


def WriteOutput(output: str, ReadDict: List[Dict[Hashable, Any]]) -> None:
    """
    Write the reads to the output file.

    Parameters
    ----------
    output : str
        The name of the output file.
    ReadDict : List[Dict[Hashable, Any]]
        A list of dictionaries, where each dictionary is a read.

    Returns
    -------
    None

    Notes
    -----
    This function takes the output file name and the dictionary of reads as input, and writes the reads
    to the output file.

    Examples
    --------
    >>> WriteOutput("output.txt", [{"Readname": "read1", "Sequence": "ATCG", "Qualities": "20"}])
    """
    with open(output, "w") as fileout:
        for index, k in enumerate(ReadDict):
            for key in ReadDict[index]:
                if key == "Readname":
                    fileout.write("@" + ReadDict[index][key] + "\n")
                if key == "Sequence":
                    fileout.write(str(ReadDict[index][key]) + "\n" + "+" + "\n")
                if key == "Qualities":
                    fileout.write(str(ReadDict[index][key]) + "\n")
