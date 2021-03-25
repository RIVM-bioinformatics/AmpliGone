import sys
import pathlib
import pandas as pd
import gzip
import pysam
import numpy as np
from io import StringIO
from Bio import SeqIO


def is_zipped(filename):
    if ".gz" in pathlib.Path(filename).suffixes:
        return True
    else:
        return False
    
def is_fastq(filename):
    fqextenstions = [".fastq", ".fq"]
    if any(item in fqextenstions for item in pathlib.Path(filename).suffixes) is True:
        return True
    else:
        return False
    
def is_bam(filename):
    if ".bam" in pathlib.Path(filename).suffixes:
        return True
    else:
        return False

def read_gzip(filename):
    with gzip.open(filename, 'rb') as f:
        lines = f.read()
    return lines

def read_fastq(filename):
    with open(filename, 'r') as f:
        lines = f.read()
    return lines

def IO_handle(inputfile):
    if is_fastq(inputfile) is True:
        ## use FastQ_IO handle
        if is_zipped(inputfile) is True:
            FastQ_IO = StringIO(read_gzip(inputfile).decode())
        else:
            FastQ_IO = StringIO(read_fastq(inputfile))
        return FastQ_IO
    elif is_bam(inputfile) is True:
        bamfile = pysam.AlignmentFile(inputfile, "rb")
        return bamfile
    else:
        print(f"\"{inputfile}\" is an unsupported filetype. Please try again with a supported filetype")
        sys.exit(-1)
        
def FlipStrand(seq, qual):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    seq = "".join(bases)
    seq = seq[::-1]
    qual = qual[::-1]
    
    return seq, qual

def ReadsToDict(inputfile):
    ReadDict = {}
    i = 0

    if is_fastq(inputfile) is True:
        handle = IO_handle(inputfile)
        for record in SeqIO.parse(handle, "fastq"):
            
            RecordQualities = "".join(
                map(lambda x: chr(x + 33), record.letter_annotations["phred_quality"])
                )
            seq = record.seq
            if len(seq) == 0:
                seq = np.nan

            ReadDict[i] = {
                "Readname": str(record.id),
                "Sequence": str(seq),
                "Qualities": str(RecordQualities)
                }
            
            i = i + 1
        handle.close()
    elif is_bam(inputfile) is True:
        handle = IO_handle(inputfile)
        for read in handle:
            if read.is_unmapped is True:
                continue
            
            readname = read.query_name
            ReadSeq = read.query_sequence
            RecordQualities = ''.join(map(lambda x: chr( x+33 ), read.query_qualities))
            
            if read.is_reverse is True:
                ReadSeq, RecordQualities = FlipStrand(ReadSeq, RecordQualities)

            if len(ReadSeq) == 0:
                ReadSeq = np.nan
            
            ReadDict[i] = {
                "Readname": str(readname), 
                "Sequence": str(ReadSeq), 
                "Qualities": str(RecordQualities)
                }
            
            i = i + 1
    
    else:
        print(f"\"{inputfile}\" is an unsupported filetype. Please try again with a supported filetype")
        sys.exit(-1)
    return ReadDict


def ReadsToList(inputfile):
    ReadList = []
    QualList = []

    if is_fastq(inputfile) is True:
        handle = IO_handle(inputfile)
        for record in SeqIO.parse(handle, "fastq"):
            ReadList.append(str(record.seq))
            quals = record.letter_annotations["phred_quality"]
            for q in quals:
                QualList.append(q)
        handle.close()
    elif is_bam(inputfile) is True:
        handle = IO_handle(inputfile)
        for read in handle:
            if read.is_unmapped is True:
                continue
            ReadSeq = read.query_sequence
            quals = read.query_qualities.tolist()
            for q in quals:
                QualList.append(q)
                
                
            ReadList.append(str(ReadSeq))
    else:
        print(f"\"{inputfile}\" is an unsupported filetype. Please try again with a supported filetype")
        sys.exit(-1)
    return ReadList, QualList

def IndexReads(inputfile):
    ReadIndex = pd.DataFrame.from_dict(ReadsToDict(inputfile), "index")
    return ReadIndex


def WriteOutput(output, ReadDict):
    with open(output, 'w') as fileout:
        for index in range(len(ReadDict)):
            for key in ReadDict[index]:
                if key == "Readname":
                    fileout.write(
                        "@" + ReadDict[index][key] + "\n"
                    )
                if key == 'Sequence':
                    fileout.write(
                        str(ReadDict[index][key]) + "\n" + "+" + "\n"
                    )
                if key == 'Qualities':
                    fileout.write(
                        str(ReadDict[index][key]) + "\n"
                    )
    pass