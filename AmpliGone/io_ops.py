import gzip
import os
import pathlib
import sys
from typing import Any, Dict, Hashable, List, TextIO, Tuple

import pandas as pd
import pgzip
import pysam
from pgzip import PgzipFile

from AmpliGone.log import log


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
    )
    primer_df = primer_df[
        ~(
            primer_df.ref.str.startswith("browser ")
            | primer_df.ref.str.startswith("track ")
        )
    ]

    return primer_df


class SequenceReads:
    def __init__(self, inputfile: str):
        log.debug(f"Starting INDEXREADS process\t@ ProcessID {os.getpid()}")
        self.tuples = []
        if self._is_fastq(inputfile):
            log.debug("INDEXREADS :: Parsing reads from FASTQ file")
            self._read_fastq(inputfile)
        elif self._is_bam(inputfile):
            log.debug("INDEXREADS :: Parsing reads from BAM file")
            self._read_bam(inputfile)
        else:
            log.error(
                f'"{inputfile}" is an unsupported filetype. Please try again with a supported filetype'
            )
            sys.exit(-1)

        log.debug("INDEXREADS :: Storing copy of indexed reads in a DataFrame")
        self.frame = pd.DataFrame.from_records(
            self.tuples, columns=["Readname", "Sequence", "Qualities"]
        )

    def _read_fastq(self, inputfile: str) -> None:
        with self._fastq_opener(inputfile) as fq:
            for line in fq:
                name = line.split()[0][1:]
                seq = next(fq).strip()
                next(fq)
                qual = next(fq).strip()

                self.tuples.append((name, seq, qual))

    def _read_bam(self, inputfile: str) -> None:
        for read in self._load_bam(inputfile):
            if read.is_unmapped:
                continue

            name: str | None = read.query_name
            seq: str | None = read.query_sequence
            qual: str | None = (
                "".join([chr(x + 33) for x in read.query_qualities])
                if read.query_qualities is not None
                else None
            )
            if name is None or seq is None or qual is None:
                continue
            if read.is_reverse:
                seq, qual = self._flip_strand(seq, qual)

            if len(seq) != len(qual):
                log.debug(
                    f"Excluding read with name '{name}' from the index due to mismatched sequence and quality length"
                )
                continue

            if len(seq) == 0:
                log.debug(
                    f"Excluding read with name '{name}' from the index due to empty sequence"
                )
                continue

            self.tuples.append((name, seq, qual))

    def _is_fastq(self, filename: str) -> bool:
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
        return any((item in ext for item in pathlib.Path(filename).suffixes))

    def _is_zipped(self, filename: str) -> bool:
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
        with open(filename, "rb") as f:
            return f.read(2) == b"\x1f\x8b"

    def _is_bam(self, filename: str) -> bool:
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
        return ".bam" in pathlib.Path(filename).suffixes

    def _load_bam(self, inputfile: str) -> pysam.AlignmentFile:
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

    def _open_gzip_fastq_file(self, filename: str) -> TextIO:
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

    def _open_fastq_file(self, filename: str) -> TextIO:
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

    def _fastq_opener(self, inputfile: str) -> TextIO:
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
        if self._is_zipped(inputfile) is True:
            log.debug("INDEXREADS :: Reading a gzipped file")
            return self._open_gzip_fastq_file(inputfile)
        log.debug("INDEXREADS :: Reading a non-gzipped file.")
        return self._open_fastq_file(inputfile)

    def _flip_strand(self, seq: str, qual: str) -> Tuple[str, str]:
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


def output_file_opener(output_file: str, threads: int) -> TextIO | PgzipFile:
    if ".gz" in output_file:
        return pgzip.open(output_file, "wt", compresslevel=6, thread=threads)
    return open(output_file, "w")


def write_output(
    output: str, read_records: List[Dict[Hashable, Any]], threads: int
) -> None:
    """
    Write the reads to the output file.

    Parameters
    ----------
    output : str
        The name of the output file.
    read_records : List[Dict[Hashable, Any]]
        A list of dictionaries, where each dictionary represents a read.
    threads : int
        The number of threads to use for writing the output file.

    Returns
    -------
    None

    Notes
    -----
    This function takes the output file name, the list of read dictionaries, and the number of threads as input.
    It writes the reads to the output file.

    Examples
    --------
    >>> write_output("output.txt", [{"Readname": "read1", "Sequence": "ATCG", "Qualities": "20"}], 4)
    """
    with output_file_opener(output, threads) as fileout:
        for index, k in enumerate(read_records):
            for key in read_records[index]:
                if key == "Readname":
                    fileout.write("@" + read_records[index][key] + "\n")
                elif key == "Sequence":
                    fileout.write(read_records[index][key] + "\n" + "+" + "\n")
                elif key == "Qualities":
                    fileout.write(read_records[index][key] + "\n")
