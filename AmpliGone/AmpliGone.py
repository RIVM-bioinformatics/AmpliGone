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
from rich import print
from rich.console import Console
from rich.progress import Progress, SpinnerColumn

from .cut_reads import CutReads
from .fasta2bed import CoordinateListsToBed, MakeCoordinateLists
from .func import FlexibleArgFormatter, RichParser, log
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

    def check_file_extensions(allowed_extensions, fname):
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
            The absolute path of the file passed.

        """

        ext = "".join(pathlib.Path(fname).suffixes)
        if not any(ext.endswith(c) for c in allowed_extensions):
            parser.error(f"File {fname} doesn't end with one of {allowed_extensions}")
        return os.path.abspath(fname)

    def check_file_exists(fname):
        """Errors if the given file `fname` does not exist

        Parameters
        ----------
        fname
            The name of the file to be read.

        Returns
        -------
            The name of the file to be read.

        """
        if os.path.isfile(fname):
            return fname
        parser.error(f'"{fname}" is not a file. Exiting...')

    parser = RichParser(
        prog="[bold]AmpliGone[/bold]",
        usage="%(prog)s \[required options] \[optional arguments]",
        description="[bold underline]AmpliGone[/bold underline]: An accurate and efficient tool to remove primers from NGS reads in reference-based experiments",
        formatter_class=FlexibleArgFormatter,
        add_help=False,
    )

    standard_threads = min(multiprocessing.cpu_count(), 128)

    required_args = parser.add_argument_group("Required Arguments")

    required_args.add_argument(
        "--input",
        "-i",
        type=lambda s: check_file_exists(
            check_file_extensions((".fastq", ".fq", ".bam", ".fastq.gz", ".fq.gz"), s)
        ),
        metavar="File",
        help="Input file with reads in either FastQ or BAM format.",
        required=True,
    )

    required_args.add_argument(
        "--output",
        "-o",
        type=lambda s: check_file_extensions((".fastq", ".fq"), s),
        metavar="File",
        help="Output (FastQ) file with cleaned reads.",
        required=True,
    )

    required_args.add_argument(
        "--reference",
        "-ref",
        type=lambda s: check_file_exists(check_file_extensions((".fasta", ".fa"), s)),
        metavar="File",
        help="Input Reference genome in FASTA format",
        required=True,
    )
    required_args.add_argument(
        "--primers",
        "-pr",
        type=lambda s: check_file_exists(
            check_file_extensions((".fasta", ".fa", ".bed"), s)
        ),
        metavar="File",
        help="""Used primer sequences in FASTA format or primer coordinates in BED format.\n Note that using bed-files overrides error-rate and ambiguity functionality""",
        required=True,
    )

    optional_args = parser.add_argument_group("Optional Arguments")

    optional_args.add_argument(
        "--amplicon-type",
        "-at",
        default="end-to-end",
        choices=("end-to-end", "end-to-mid", "fragmented"),
        help="Define the amplicon-type, either being [green]'end-to-end'[/green], [green]'end-to-mid'[/green], or [green]'fragmented'[/green]. See the docs for more info :book:",
        required=False,
        metavar="'end-to-end'/'end-to-mid'/'fragmented'",
    )

    optional_args.add_argument(
        "--fragment-lookaround-size",
        "-fls",
        required=False,
        type=int,
        metavar="N",
        help="The number of bases to look around a primer-site to consider it part of a fragment. Only used if amplicon-type is 'fragmented'. Default is 10",
    )

    optional_args.add_argument(
        "--export-primers",
        "-ep",
        type=lambda s: check_file_extensions((".bed",), s),
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
        help=f"""Number of threads you wish to use.\n Default is the number of available threads in your system ({standard_threads})""",
    )

    optional_args.add_argument(
        "-to",
        action="store_true",
        help="If set, AmpliGone will always create the output files even if there is nothing to output. (for example when an empty input-file is given)\n This is useful in (automated) pipelines where you want to make sure that the output files are always created.",
        required=False,
    )

    optional_args.add_argument(
        "--error-rate",
        "-er",
        type=float,
        default=0.1,
        metavar="N",
        help="The maximum allowed error rate for the primer search. Use 0 for exact primer matches.",
        required=False,
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
    fragment_lookaround_size,
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
            fragment_lookaround_size,
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
    log.info(
        f"Starting AmpliGone with inputfile [green]'{os.path.abspath(args.input)}'[/green]"
    )

    with cf.ThreadPoolExecutor(max_workers=args.threads) as ex:
        TP_indexreads = ex.submit(IndexReads, args.input)

        if not args.primers.endswith(".bed"):
            TP_PrimerLists = ex.submit(
                MakeCoordinateLists, args.primers, args.reference, args.error_rate
            )
            primer_df = TP_PrimerLists.result()
        else:
            log.info(
                "Primer coordinates are given in BED format, skipping primer search"
            )
            primer_df = read_bed(args.primers)
        IndexedReads = TP_indexreads.result()

    if len(IndexedReads.index) < 1 and args.to is True:
        ReadDict = IndexedReads.to_dict(orient="records")
        WriteOutput(args.output, ReadDict)
        if args.export_primers is not None:
            with open(args.export_primers, "w") as f:
                f.write("")
        log.warning(
            "AmpliGone was given an empty input file but the [green]'-to'[/green] flag was given.\nOne or multiple empty output file(s) have therefore been generated.\n[bold yellow]Please check the input file to make sure this is correct[/bold yellow]"
        )
        sys.exit(0)
    elif len(IndexedReads.index) < 1:
        log.error(
            "AmpliGone was given an empty input file. Exiting..\nPlease check the input file to make sure this is correct\n\n[bold yellow]Use the -to flag to force AmpliGone to create an output file even if there is nothing to output.[/bold yellow]"
        )
        sys.exit(1)
    else:
        log.info(
            f"Succesfully loaded [bold green]{len(IndexedReads.index)}[/bold green] reads."
        )

    if len(primer_df) < 1:
        log.error(
            "AmpliGone was unable to match any primers to the reference. AmpliGone is therefore unable to remove primers from the reads.\nPlease check the primers and reference to make sure this is correct\nExiting..."
        )
        sys.exit(1)

    if len(IndexedReads.index) < args.threads:
        log.info(
            f"[yellow]AmpliGone is set to use more threads than reads present. Downscaling threads to match.[/yellow]"
        )
        args.threads = len(IndexedReads.index)

    log.info("Finding optimal alignment parameters for the given reads")
    # Todo: split this over two threads if possible
    if len(IndexedReads.index) > 20000:
        preset, scoring = FindPreset(
            args.threads, IndexedReads.sample(frac=0.3)
        )  # Todo: Make this more efficient
    else:
        preset, scoring = FindPreset(args.threads, IndexedReads)

    ## correct the lookaround size if the amplicon type is not fragmented
    if args.amplicon_type != "fragmented":
        args.fragment_lookaround_size = 10000
    elif args.fragment_lookaround_size is None:
        args.fragment_lookaround_size = 10
        log.warning(
            "[yellow]No fragment lookaround size was given, [underline]using default of 10[/underline][/yellow]"
        )

    log.info(
        f"Distributing {len(IndexedReads.index)} reads across {args.threads} threads for processing. Processing around [bold green]{round(len(IndexedReads.index)/args.threads)}[/bold green] reads per thread"
    )

    IndexedReads = IndexedReads.sample(frac=1).reset_index(drop=True)

    with Progress(
        SpinnerColumn(),
        *Progress.get_default_columns(),
        console=Console(record=True),
        transient=True,
    ) as progress:
        progress.add_task("[yellow]Removing primer sequences...", total=None)
        ProcessedReads = parallel(
            IndexedReads,
            CutReads,
            args.threads,
            primer_df,
            args.reference,
            preset,
            scoring,
            fragment_lookaround_size=args.fragment_lookaround_size,
            amplicon_type=args.amplicon_type,
        )

        ProcessedReads.reset_index(drop=True)
    log.info("Done removing primer sequences")

    total_nuc_preprocessing = sum(tuple(chain(IndexedReads["Sequence"].str.len())))
    removed_coordinates = tuple(chain(*ProcessedReads["Removed_coordinates"]))

    log.info(
        f"\tRemoved a total of [bold cyan]{len(removed_coordinates)}[/bold cyan] nucleotides."
    )
    log.info(
        f"\tThis is [bold cyan]{round((len(removed_coordinates)/total_nuc_preprocessing)*100, 2)}%[/bold cyan] of the total amount of nucleotides present in the original reads."
    )
    log.info(
        f"\tThese nucleotides were removed from [bold cyan]{len(set(removed_coordinates))}[/bold cyan] unique nucleotide-coordinates."
    )
    log.info(
        f"\tThese nucleotide-coordinates correspond to the coordinates of [bold cyan]{len(primer_df)}[/bold cyan] (found) primers."
    )

    log.info("Writing output files")
    if args.export_primers is not None:
        removed_coords = set(removed_coordinates)
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
