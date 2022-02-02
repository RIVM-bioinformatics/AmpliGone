import gzip
import pathlib
import sys

import pandas as pd
import pysam


def is_zipped(filename):
    return bool(".gz" in pathlib.Path(filename).suffixes)


def is_fastq(filename):
    ext = [".fastq", ".fq"]
    return bool(any(item in ext for item in pathlib.Path(filename).suffixes))


def is_bam(filename):
    return bool(".bam" in pathlib.Path(filename).suffixes)


def read_gzip(filename):
    return gzip.open(filename, "rt")


def read_fastq(filename):
    return open(filename, "rt")


def fastq_opener(inputfile):
    if is_zipped(inputfile) is True:
        return read_gzip(inputfile)
    return read_fastq(inputfile)


def LoadBam(inputfile):
    return pysam.AlignmentFile(inputfile, "rb")


def FlipStrand(seq, qual):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    seq = "".join(bases)
    seq = seq[::-1]
    qual = qual[::-1]

    return seq, qual


def LoadData(inputfile):
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
    print(
        f'"{inputfile}" is an unsupported filetype. Please try again with a supported filetype'
    )
    sys.exit(-1)


def IndexReads(inputfile):
    return pd.DataFrame.from_records(
        LoadData(inputfile), columns=["Readname", "Sequence", "Qualities"]
    )


def WriteOutput(output, ReadDict):
    with open(output, "w") as fileout:
        for index, k in enumerate(ReadDict):
            for key in ReadDict[index]:
                if key == "Readname":
                    fileout.write("@" + ReadDict[index][key] + "\n")
                if key == "Sequence":
                    fileout.write(str(ReadDict[index][key]) + "\n" + "+" + "\n")
                if key == "Qualities":
                    fileout.write(str(ReadDict[index][key]) + "\n")
