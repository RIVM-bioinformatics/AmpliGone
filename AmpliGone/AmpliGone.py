"""
Copyright © 2021 RIVM

https://github.com/RIVM-bioinformatics/AmpliGone
"""


import argparse
import concurrent.futures as cf
import multiprocessing
import os
import pathlib
import sys
from itertools import chain

import numpy as np
import pandas as pd
import parmap

from .cut_reads import CutReads
from .fasta2bed import CoordinateListsToBed, MakeCoordinateLists
from .func import MyHelpFormatter, color
from .io_ops import IndexReads, WriteOutput, read_bed
from .mappreset import FindPreset
from .version import __version__


def get_args(givenargs):
    """It takes the arguments given to the script and parses them into the argparse namespace

    Parameters
    ----------
    givenargs
        the arguments given to the script

    Returns
    -------
        The arguments that are given to the program.

    """

    def fastq_or_bam(choices, fname):
        """If the input file exists, check that it ends with one of the extensions in the list of choices. If
        it does, return the file name. If it doesn't, print an error message and exit

        Parameters
        ----------
        choices
            a list of file extensions that are allowed
        fname
            the name of the file to be processed

        Returns
        -------
            the file name if it is a file and if it ends with one of the choices.

        """
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffixes)
            if ext not in choices:
                parser.error(f"Input file doesn't end with one of {choices}")
            return fname
        print(f'"{fname}" is not a file. Exiting...')
        sys.exit(-1)

    def fastq_output(choices, fname):
        """If the file extension of the input file is not one of the choices, then raise an error

        Parameters
        ----------
        choices
            a list of file extensions that are allowed
        fname
            The name of the file to be read.

        Returns
        -------
            The file name.

        """
        ext = "".join(pathlib.Path(fname).suffixes)
        if ext not in choices:
            parser.error(f"Input file doesn't end with one of {choices}")
        return fname

    def check_extensions(choices, fname):
        """It checks that the file extension of the file name passed to it is one of the extensions in the list
        passed to it

        Parameters
        ----------
        choices
            a list of strings that are valid extensions
        fname
            The name of the file to be read.

        Returns
        -------
            The file name.

        """
        ext = "".join(pathlib.Path(fname).suffixes)
        if ext not in choices:
            parser.error(
                f"File doesn't end with {choices[0] if len(choices) == 1 else f'one of {choices}'}"
            )
        return fname

    parser = argparse.ArgumentParser(
        prog="AmpliGone",
        usage="%(prog)s [required options] [optional arguments]",
        description="AmpliGone: An accurate and efficient tool to remove primers from NGS reads in reference-based experiments",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )

    standard_threads = min(multiprocessing.cpu_count(), 128)

    required_args = parser.add_argument_group("Required Arguments")

    required_args.add_argument(
        "--input",
        "-i",
        type=lambda s: fastq_or_bam(
            (".fastq", ".fq", ".bam", ".fastq.gz", ".fq.gz"), s
        ),
        metavar="File",
        help="Input file with reads in either FastQ or BAM format.",
        required=True,
    )

    required_args.add_argument(
        "--output",
        "-o",
        type=lambda s: fastq_output((".fastq", ".fq"), s),
        metavar="File",
        help="Output (FastQ) file with cleaned reads.",
        required=True,
    )

    required_args.add_argument(
        "--reference",
        "-ref",
        type=lambda s: check_extensions((".fasta", ".fa"), s),
        metavar="File",
        help="Input Reference genome in FASTA format",
        required=True,
    )
    required_args.add_argument(
        "--primers",
        "-pr",
        type=lambda s: check_extensions((".fasta", ".fa", ".bed"), s),
        metavar="File",
        help="Used primer sequences in FASTA format or primer coordinates in BED format.\nNote that using bed-files overrides error-rate and ambiguity functionality",
        required=True,
    )

    optional_args = parser.add_argument_group("Optional Arguments")
    
    optional_args.add_argument(
        "--amplicon-type",
        "-at",
        default="end-to-end",
        choices=("end-to-end", "end-to-mid"),
        help="Define the amplicon-type, either being 'end-to-end' or 'end-to-mid'.\nSee the docs for more info",
        required=False,
    )

    optional_args.add_argument(
        "--export-primers",
        "-ep",
        type=lambda s: check_extensions((".bed",), s),
        metavar="File",
        help="Output BED file with found primer coordinates if they are actually cut from the reads",
        required=False,
    )

    optional_args.add_argument(
        "--threads",
        "-t",
        type=int,
        default=standard_threads,
        metavar="N",
        help=f"Number of threads you wish to use.\nDefault is the number of available threads in your system ({standard_threads})",
    )

    optional_args.add_argument(
        "--version",
        "-v",
        action="version",
        version=__version__,
        help="Show the AmpliGone version and exit",
    )

    optional_args.add_argument(
        "--help",
        "-h",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )

    optional_args.add_argument(
        "-to",
        action="store_true",
        help="If set, AmpliGone will always create the output files even if there is nothing to output. (for example when an empty input-file is given)\nThis is useful in (automated) pipelines where you want to make sure that the output files are always created.",
        required=False,
    )

    optional_args.add_argument(
        "--error-rate",
        "-er",
        type=float,
        default=0.1,
        metavar="N",
        help="The maximum allowed error rate for the primer search. Use 0 for exact primer matches.\nDefault is '0.1'. ",
        required=False,
    )

    flags = parser.parse_args(givenargs)

    return flags


def parallel(
    frame,
    function,
    workers,
    primer_df,
    reference,
    preset,
    scoring,
    amplicon_type,
):
    frame_split = np.array_split(frame, workers)
    tr = [*range(workers)]
    readframe = pd.concat(
        parmap.map(
            function,
            zip(frame_split, tr),
            primer_df,
            reference,
            preset,
            scoring,
            amplicon_type,
            workers,
            pm_processes=workers,
        )
    )
    return readframe


def main():
    if len(sys.argv[1:]) < 1:
        print(
            "AmpliGone was called but no arguments were given, please try again.\nUse 'AmpliGone -h' to see the help document"
        )
        sys.exit(1)
    args = get_args(sys.argv[1:])

    with cf.ThreadPoolExecutor(max_workers=args.threads) as ex:
        TP_indexreads = ex.submit(IndexReads, args.input)

        if not args.primers.endswith(".bed"):
            TP_PrimerLists = ex.submit(
                MakeCoordinateLists, args.primers, args.reference, args.error_rate
            )
            primer_df = TP_PrimerLists.result()
        else:
            primer_df = read_bed(args.primers)
        # print(primer_df)
        # exit(0)
        IndexedReads = TP_indexreads.result()

    if len(IndexedReads.index) < 1 and args.to is True:
        ReadDict = IndexedReads.to_dict(orient="records")
        WriteOutput(args.output, ReadDict)
        if args.export_primers is not None:
            with open(args.export_primers, "w") as f:
                f.write("")
        print(
            f"""
    {color.RED}AmpliGone was given an empty input file but the '-to' flag was given.
    {color.YELLOW}One or multiple empty output file(s) have therefore been generated.
    {color.RED}Please check the input file to make sure this is correct{color.END}
    """
        )
        sys.exit(0)
    elif len(IndexedReads.index) < 1:
        print(
            f"""
    {color.RED}AmpliGone was given an empty input file. Exiting...
    Please check the input file to make sure this is correct{color.END}

    {color.YELLOW}Use the -to flag to force AmpliGone to create an output file even if there is nothing to output.{color.END}
    """
        )
        sys.exit(1)
    else:
        print(
            f"""
    Succesfully loaded {color.BOLD + color.GREEN}{len(IndexedReads.index)}{color.END} reads.
            """
        )

    if len(primer_df) < 1:
        print(
            f"""
    {color.RED}AmpliGone was unable to match any primers to the reference. AmpliGone is therefore unable to remove primers from the reads.
    {color.RED}Please check the primers and reference to make sure this is correct
    {color.RED}Exiting...{color.END}
    """
        )
        sys.exit(1)

    # Todo: split this over two threads if possible
    if len(IndexedReads.index) > 20000:
        preset, scoring = FindPreset(
            args.threads, IndexedReads.sample(frac=0.3)
        )  # Todo: Make this more efficient
    else:
        preset, scoring = FindPreset(args.threads, IndexedReads)

    IndexedReads = IndexedReads.sample(frac=1).reset_index(drop=True)

    ProcessedReads = parallel(
        IndexedReads,
        CutReads,
        args.threads,
        primer_df,
        args.reference,
        preset,
        scoring,
        amplicon_type=args.amplicon_type,
    )
    ProcessedReads.reset_index(drop=True)

    if args.export_primers is not None:
        removed_coords = set(chain(*ProcessedReads["Removed_coordinates"]))
        filtered_primer_df = primer_df[
            primer_df[["start", "end"]].apply(
                lambda r: any(coord in removed_coords for coord in range(*r)),
                axis=1,
            )
        ]
        CoordinateListsToBed(filtered_primer_df, args.export_primers)

    ProcessedReads = ProcessedReads.drop(columns=["Removed_coordinates"])

    ReadDict = ProcessedReads.to_dict(orient="records")

    WriteOutput(args.output, ReadDict)
