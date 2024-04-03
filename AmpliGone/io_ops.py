import gzip
import pathlib
import sys

import pandas as pd
import pysam

from .func import log


def is_zipped(filename):
    """If the filename has a ".gz" suffix, return True, otherwise return False

    Parameters
    ----------
    filename
        the name of the file to be read

    Returns
    -------
        A boolean value.

    """
    return bool(".gz" in pathlib.Path(filename).suffixes)


def is_fastq(filename):
    """If any of the items in the list `ext` are in the list of suffixes of the file `filename`, then
    return `True`, otherwise return `False`

    Parameters
    ----------
    filename
        The name of the file to be checked.

    Returns
    -------
        A boolean value.

    """
    ext = [".fastq", ".fq"]
    return bool(any(item in ext for item in pathlib.Path(filename).suffixes))


def is_bam(filename):
    """It returns True if the filename ends in ".bam", return False otherwise

    Parameters
    ----------
    filename
        the name of the file to be checked

    Returns
    -------
        A boolean value.

    """
    return bool(".bam" in pathlib.Path(filename).suffixes)


def read_gzip(filename):
    """It opens a gzip file for reading, and returns an opened file object.

    Parameters
    ----------
    filename
        The name of the file to open.

    Returns
    -------
        A file object

    """
    return gzip.open(filename, "rt")


def read_fastq(filename):
    """It opens a file and returns a file object

    Parameters
    ----------
    filename
        the name of the file to read

    Returns
    -------
        An opened file object

    """
    return open(filename, "rt")


def fastq_opener(inputfile):
    """It opens a fastq file, and if it's zipped, it unzips it

    Parameters
    ----------
    inputfile
        the name of the input file

    Returns
    -------
        A file object

    """
    if is_zipped(inputfile) is True:
        return read_gzip(inputfile)
    return read_fastq(inputfile)


def LoadBam(inputfile):
    """It loads a BAM file and returns a pysam.AlignmentFile object

    Parameters
    ----------
    inputfile
        the path to the BAM file

    Returns
    -------
        A pysam.AlignmentFile object

    """
    return pysam.AlignmentFile(inputfile, "rb")


def read_bed(filename):
    """It reads a BED file, removes the original BEDfile header lines, and returns a pandas dataframe

    Parameters
    ----------
    filename
        the name of the file needs to be read

    Returns
    -------
        A dataframe with the columns ref, start, end, name, score, and strand.

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


def FlipStrand(seq, qual):
    """It takes a sequence and a quality score, and returns the reverse complement of the sequence and the
    reverse of the quality score

    Parameters
    ----------
    seq
        the sequence of the read
    qual
        the quality score of the read

    Returns
    -------
        the sequence and quality score in reverse complement order.

    """
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    seq = "".join(bases)
    seq = seq[::-1]
    qual = qual[::-1]

    return seq, qual


def LoadData(inputfile):
    """It takes a file and returns a list of tuples, where each tuple contains the name, sequence, and
    quality score of a read

    Parameters
    ----------
    inputfile
        The input file to be processed.

    Returns
    -------
        A list of tuples. Each tuple contains the name, sequence, and quality score of a read.

    """
    Reads = []

    if is_fastq(inputfile) is True:
        with fastq_opener(inputfile) as fq:
            for line in fq:
                name = line.split()[0][1:]
                seq = next(fq).strip()
                next(fq)
                qual = next(fq).strip()

                Reads.append((name, seq, qual))
        return Reads
    if is_bam(inputfile) is True:
        for read in LoadBam(inputfile):
            if read.is_unmapped is True:
                continue

            name = read.query_name
            seq = read.query_sequence
            qual = "".join(map(lambda x: chr(x + 33), read.query_qualities))

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


def IndexReads(inputfile):
    """It reads in an input file and returns a dataframe with the readname, sequence, and qualities

    Parameters
    ----------
    inputfile
        The path to the input file.

    Returns
    -------
        A dataframe with the columns Readname, Sequence, and Qualities.

    """
    return pd.DataFrame.from_records(
        LoadData(inputfile), columns=["Readname", "Sequence", "Qualities"]
    )


def WriteOutput(output, ReadDict):
    """This function takes the output file name and the dictionary of reads as input, and writes the reads
    to the output file

    Parameters
    ----------
    output
        the name of the output file
    ReadDict
        a dictionary of dictionaries, where each dictionary is a read.

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
