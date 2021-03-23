"""
>> insert copyright notices here
"""


import os
import sys
import argparse
import pathlib
import multiprocessing
import time
import parmap
import numpy as np
import pandas as pd
import concurrent.futures as cf

from .version import __version__
from .func import MyHelpFormatter, color
from .io_ops import IndexReads, WriteOutput
from .CoordinateSearch import MakeCoordinateLists
from .mappreset import FindPreset
from .cut_reads import End_to_End, End_to_Mid
    
    
def get_args(givenargs):
    """
    Parse the cmd args
    """
    
    def fastq_or_bam(choices,fname):
        if os.path.isfile(fname):
            ext = ''.join(pathlib.Path(fname).suffixes)
            if ext not in choices:
                parser.error("Input file doesn't end with one of {}".format(choices))
            return fname
        else:
            print(f'\"{fname}\" is not a file. Exiting...')
            sys.exit(-1)
    
    def fastq_output(choices,fname):
        ext = ''.join(pathlib.Path(fname).suffixes)
        if ext not in choices:
            parser.error("Input file doesn't end with one of {}".format(choices))
        return fname
    
    def fasta_input(choices,fname):
        if os.path.isfile(fname):
            ext = ''.join(pathlib.Path(fname).suffixes)
            if ext not in choices:
                parser.error("Input file doesn't end with one of {}".format(choices))
            return fname
        else:
            print(f'\"{fname}\" is not a file. Exiting...')
            sys.exit(-1)
    
    parser = argparse.ArgumentParser(
        prog='AmpliGone',
        usage='%(prog)s [required options] [optional arguments]',
        description="AmpliGone: An accurate and efficient tool to remove primers from NGS reads in reference-based experiments",
        formatter_class=MyHelpFormatter, add_help=False
    )
    
    standard_threads = min(multiprocessing.cpu_count(), 128)
    
    parser.add_argument(
        '--input', '-i',
        type=lambda s:fastq_or_bam((".fastq", ".fq", ".bam", ".fastq.gz", ".fq.gz"),s),
        metavar='File',
        help="Input file with reads in either FastQ or BAM format.",
        required=True
    )
    
    parser.add_argument(
        '--output', '-o',
        type=lambda s:fastq_output((".fastq", ".fq"),s),
        metavar='File',
        help='Output (FastQ) file with cleaned reads.',
        required=True
    )
    
    parser.add_argument(
        '--reference', '-ref',
        type=lambda s:fasta_input((".fasta", ".fa"),s),
        metavar='File',
        help='Input Reference genome in FASTA format',
        required=True
    )
    parser.add_argument(
        '--primers', '-pr',
        type=lambda s:fasta_input((".fasta", ".fa"),s),
        metavar='File',
        help="Used primer sequences in FASTA format",
        required=True
    )
    parser.add_argument(
        '--amplicon-type',
        '-at',
        const='end-to-end',
        nargs='?',
        choices=('end-to-end', 'end-to-mid'),
        help="Define the amplicon-type, either being 'end-to-end' or 'end-to-mid'.\nSee the docs for more info",
        required=True
    )
    parser.add_argument(
        '--threads', '-t',
        type=int,
        default=standard_threads,
        metavar="N",
        help=f"Number of threads you wish to use.\nDefault is the number of available threads in your system ({standard_threads})"
    )
    
    parser.add_argument(
        '--version', '-v',
        action='version',
        version=__version__,
        help='Show the AmpliGone version and exit'
    )
    
    parser.add_argument(
        '--help', '-h',
        action='help',
        default=argparse.SUPPRESS,
        help="Show this help message and exit"
    )
    
    flags = parser.parse_args(givenargs)
    
    return flags

def parallel(frame, function, workers, LeftPrimers, RightPrimers, reference, preset):
    frame_split = np.array_split(frame, workers)
    tr = [*range(workers)]
    readframe = pd.concat(parmap.map(function, zip(frame_split, tr), LeftPrimers, RightPrimers, reference, preset, workers, pm_processes=workers))
    return readframe

def main():
    if len(sys.argv[1:]) < 1:
        print('AmpliGone was called but no arguments were given, please try again.\nUse \'AmpliGone -h\' to see the help document')
        sys.exit(1)
    args = get_args(sys.argv[1:])
    
    with cf.ThreadPoolExecutor(max_workers=args.threads) as exec:
        TP_indexreads = exec.submit(IndexReads, args.input)
        TP_PrimerLists = exec.submit(MakeCoordinateLists, args.primers, args.reference)
        TP_FindPreset = exec.submit(FindPreset, args.input, args.threads)
        
        IndexedReads = TP_indexreads.result()
        LeftPrimers, RightPrimers, CombPrimers = TP_PrimerLists.result()
        preset = TP_FindPreset.result()
    
    IndexedReads.dropna(subset=['Sequence'], inplace=True)
    IndexedReads = IndexedReads.sample(frac=1).reset_index(drop=True)
    
    if args.amplicon_type == "end-to-end":
        
        ProcessedReads = parallel(IndexedReads, End_to_End, args.threads, LeftPrimers, RightPrimers, args.reference, preset)
        ProcessedReads.reset_index(drop=True)

    if args.amplicon_type == "end-to-mid":
        ProcessedReads = parallel(IndexedReads, End_to_Mid, args.threads, LeftPrimers, RightPrimers, args.reference, preset)
        ProcessedReads.reset_index(drop=True)
        
    ReadDict = ProcessedReads.to_dict(orient='records')
    
    WriteOutput(args.output, ReadDict)