"""
Copyright © 2021 RIVM

https://github.com/RIVM-bioinformatics/AmpliGone
"""

import argparse
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import chain
from typing import Callable

import pandas as pd
import parmap
from rich import print as pprint
from rich.console import Console
from rich.progress import Progress, SpinnerColumn

import AmpliGone.alignmentmatrix as AlignmentMatrix
import AmpliGone.alignmentpreset as AlignmentPreset
from AmpliGone import __prog__, __version__
from AmpliGone.args import get_args
from AmpliGone.cut_reads import cut_reads
from AmpliGone.fasta2bed import coord_lists_to_bed, find_or_read_primers
from AmpliGone.io_ops import SequenceReads, write_output
from AmpliGone.log import log


def check_loaded_index(
    indexed_reads: SequenceReads, args: argparse.Namespace
) -> SequenceReads:
    """
    Check the input file for indexed reads and handle empty input file scenarios.

    Parameters
    ----------
    indexed_reads : SequenceReads
        The indexed reads.

    args : argparse.Namespace
        The command-line arguments.

    Returns
    -------
    SequenceReads
        The indexed reads.

    Raises
    ------
    SystemExit
        If the input file is empty and the '-to' flag is not provided.

    Notes
    -----
    This function checks if the input file for indexed reads is empty. If it is empty and the '-to' flag is not provided, it raises a SystemExit exception. If the '-to' flag is provided, it writes an empty output file and logs a warning message. If the input file is not empty, it logs an info message and returns the indexed reads.

    """
    if len(indexed_reads.tuples) < 1:
        if args.to is True:
            read_records = indexed_reads.frame.to_dict(orient="records")
            write_output(args.output, read_records, threads=args.threads)
            if args.export_primers is not None:
                with open(args.export_primers, "w", encoding="utf-8") as f:
                    f.write("")
            log.warning(
                f"{__prog__} was given an empty input file but the [green]'-to'[/green] flag was given.\nOne or multiple empty output file(s) have therefore been generated.\n[bold yellow]Please check the input file to make sure this is correct[/bold yellow]"
            )
            sys.exit(0)
        else:
            log.error(
                f"{__prog__} was given an empty input file. Exiting..\nPlease check the input file to make sure this is correct\n\n[bold yellow]Use the -to flag to force {__prog__} to create an output file even if there is nothing to output.[/bold yellow]"
            )
            sys.exit(1)
    log.info(
        f"Succesfully loaded [bold green]{len(indexed_reads.tuples)}[/bold green] reads."
    )
    return indexed_reads


def primer_df_to_primer_index(
    primer_df: pd.DataFrame, bind_virtual_primer: bool = True
) -> tuple[defaultdict, defaultdict]:
    """
    Convert primer DataFrame to primer index dictionaries.

    Parameters
    ----------
    primer_df : pd.DataFrame
        DataFrame containing primer information.
    bind_virtual_primer : bool, optional
        Whether to match closely positioned primers in the same orientation into a virtual primer. Defaults to False.

    Returns
    -------
    tuple[defaultdict, defaultdict]
        A tuple of two defaultdicts representing the forward and reverse primer indices.

    Raises
    ------
    SystemExit
        If the primer DataFrame is empty.

    Notes
    -----
    This function converts the primer DataFrame into two defaultdicts representing the forward and reverse primer indices. The primer DataFrame should contain information about the primers, including the reference ID, start coordinate, end coordinate, and strand. If the primer DataFrame is empty, a SystemExit exception is raised.

    The bind_virtual_primer parameter determines whether virtual primers should be bound. If bind_virtual_primer is set to True, closely positioned primers in the same orientation will be matched to generate a virtual primer for proper removal. By default, bind_virtual_primer is set to False.

    Examples
    --------
    >>> primer_df = pd.DataFrame(...)
    >>> forward_dict, reverse_dict = primer_df_to_primer_index(primer_df)
    >>> print(forward_dict)
    defaultdict(<class 'set'>, {'ref1': {1, 2, 3}, 'ref2': {4, 5, 6}})
    >>> print(reverse_dict)
    defaultdict(<class 'set'>, {'ref1': {7, 8, 9}, 'ref2': {10, 11, 12}})
    """
    if len(primer_df) < 1:
        log.error(
            f"{__prog__} was unable to match any primers to the reference. {__prog__} is therefore unable to remove primers from the reads.\nPlease check the primers and reference to make sure this is correct\nExiting..."
        )
        sys.exit(1)

    forward_dict: defaultdict[str, set[None | int]] = defaultdict(set)
    reverse_dict: defaultdict[str, set[None | int]] = defaultdict(set)

    reference_set: set[str] = set(primer_df["ref"].unique())
    for refid in reference_set:
        forward_dict[refid] = set()
        reverse_dict[refid] = set()

    # split the primer_df into forward and reverse primers based on strand
    forward_primers_df, reverse_primer_df = (
        primer_df[primer_df["strand"] == "+"],
        primer_df[primer_df["strand"] == "-"],
    )

    def coordinates_to_index(
        oriented_primer_df: pd.DataFrame, refid: str, start: int, end: int
    ) -> range:
        """
        Converts coordinates to an index range based on the given parameters.

        Parameters
        ----------
        oriented_primer_df : pd.DataFrame
            The DataFrame containing the oriented primer information.
        refid : str
            The reference ID.
        start : int
            The start coordinate.
        end : int
            The end coordinate.

        Returns
        -------
        range
            The index range based on the given coordinates.

        Notes
        -----
        If bind_virtual_primer is True, the function will also match closely positioned primers in the same orientation to generate a virtual primer for proper removal.
        bind_virtual_primer is set in the wrapping function and inherited here.
        """
        if bind_virtual_primer:
            length = end - start
            for _, refid2, start2, end2 in oriented_primer_df[
                ["ref", "start", "end"]
            ].itertuples():
                if (
                    refid == refid2
                    and start2 <= end + length
                    and end2 >= start - length
                ):
                    start = min(start, start2)
                    end = max(end, end2)
        return range(start + 1, end)

    # iterate over forward_primers_df and add the coordinates between "start" and "end" to the forward_dict
    for _, refid, start, end in forward_primers_df[
        ["ref", "start", "end"]
    ].itertuples():
        forward_dict[refid].update(
            coordinates_to_index(forward_primers_df, refid, start, end)
        )

    # iterate over reverse_primers_df and add the coordinates between "start" and "end" to the reverse_dict
    for _, refid, start, end in reverse_primer_df[["ref", "start", "end"]].itertuples():
        reverse_dict[refid].update(
            coordinates_to_index(reverse_primer_df, refid, start, end)
        )

    return forward_dict, reverse_dict


def check_thread_count(
    indexed_reads: SequenceReads, args: argparse.Namespace
) -> argparse.Namespace:
    """
    Check the thread count and adjust it if necessary based on the number of reads.

    Parameters
    ----------
    indexed_reads : SequenceReads
        The indexed reads.
    args : argparse.Namespace
        The command-line arguments.

    Returns
    -------
    argparse.Namespace
        The updated command-line arguments.

    Notes
    -----
    This function checks if the number of threads specified in the command-line arguments
    is greater than the number of reads present in the indexed reads. If it is, the number
    of threads is downscaled to match the number of reads.

    If the number of reads is less than the specified number of threads, the number of threads
    is set to the minimum value between the number of reads and 2.

    Examples
    --------
    >>> indexed_reads = SequenceReads(...)
    >>> args = argparse.Namespace(threads=4)
    >>> updated_args = check_thread_count(indexed_reads, args)
    >>> updated_args.threads
    2
    """
    if len(indexed_reads.tuples) < args.threads:
        log.info(
            f"[yellow]{__prog__} is set to use more threads than reads present. Downscaling threads to match.[/yellow]"
        )
        args.threads = min(len(indexed_reads.tuples), 2)
    return args


def correct_fragment_lookaround_size(args: argparse.Namespace) -> argparse.Namespace:
    """
    Corrects the fragment lookaround size based on the amplicon type.

    Parameters
    ----------
    args : argparse.Namespace
        The command-line arguments.

    Returns
    -------
    argparse.Namespace
        The updated command-line arguments.

    Notes
    -----
    This function adjusts the `fragment_lookaround_size` based on the `amplicon_type` value in the `args` namespace.
    If the `amplicon_type` is not "fragmented", the `fragment_lookaround_size` is set to 10000.
    If the `amplicon_type` is "fragmented" and `fragment_lookaround_size` is not provided, it is set to 10.
    A warning message is logged if the `fragment_lookaround_size` is set to the default value.

    Examples
    --------
    >>> args = argparse.Namespace(amplicon_type="fragmented", fragment_lookaround_size=None)
    >>> corrected_args = correct_fragment_lookaround_size(args)
    >>> print(corrected_args.fragment_lookaround_size)
    10
    """
    if args.amplicon_type != "fragmented":
        args.fragment_lookaround_size = 10000
    elif args.fragment_lookaround_size is None:
        args.fragment_lookaround_size = 10
        log.warning(
            "[yellow]No fragment lookaround size was given, [underline]using default of 10[/underline][/yellow]"
        )
    return args


def parallel_dispatcher(
    indexed_reads: SequenceReads,
    args: argparse.Namespace,
    primer_sets: tuple[defaultdict, defaultdict],
    preset: str,
    matrix: list[int],
) -> pd.DataFrame:
    """
    Wrapping function that actually calls the parallelization function to process the primer removal process of the reads.

    Parameters
    ----------
    indexed_reads : SequenceReads
        The indexed reads to be processed.
    args : argparse.Namespace
        The command-line arguments.
    primer_sets : tuple[defaultdict, defaultdict]
        The primer sequences to be removed.
    preset : str
        The preset configuration for processing.
    matrix : list[int]
        The matrix for processing.

    Returns
    -------
    pd.DataFrame
        The processed reads with primer sequences removed.
    """
    with Progress(
        SpinnerColumn(),
        *Progress.get_default_columns(),
        console=Console(
            record=False,
        ),
        transient=True,
        disable=args.quiet or args.verbose,
    ) as progress:
        progress.add_task("[yellow]Removing primer sequences...", total=None)
        processed_reads = parallel(
            indexed_reads.frame,
            cut_reads,
            args.threads,
            primer_sets,
            args.reference,
            preset,
            matrix,
            fragment_lookaround_size=args.fragment_lookaround_size,
            amplicon_type=args.amplicon_type,
        )

        processed_reads.reset_index(drop=True)
    log.info("Done removing primer sequences")
    return processed_reads


def parallel(
    frame: pd.DataFrame,
    function: Callable[..., pd.DataFrame],
    workers: int,
    primer_sets: tuple[defaultdict, defaultdict],
    reference: str,
    preset: str,
    scoring: list[int],
    fragment_lookaround_size: int,
    amplicon_type: str,
) -> pd.DataFrame:
    """
    Apply a function to a pandas DataFrame in parallel using multiple workers.

    Parameters
    ----------
    frame : pandas.DataFrame
        The DataFrame to apply the function to.
    function : Callable[..., pandas.DataFrame]
        The function to apply to the DataFrame.
    workers : int
        The number of workers to use for parallel processing.
    primer_df : tuple[defaultdict, defaultdict]
        A tuple containing the indexes of the primer coordinates to remove.
    reference : str
        The reference sequence to use for alignment.
    preset : str
        The preset to use for alignment.
    scoring : list[int]
        The scoring matrix to use for alignment.
        The size of the fragment lookaround.
    fragment_lookaround_size : int
    amplicon_type : str
        The type of amplicon.

    Returns
    -------
    pandas.DataFrame
        The resulting DataFrame after applying the function to the input DataFrame in parallel.
    """
    frame_split = [frame.iloc[i::workers] for i in range(workers)]
    tr = [*range(workers)]
    df = pd.concat(
        parmap.map(
            function,
            zip(frame_split, tr),
            primer_sets,
            reference,
            preset,
            scoring,
            fragment_lookaround_size,
            amplicon_type,
            pm_processes=workers,
        )
    )
    # parmap.map sometimes returns Any, but we know it's a DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"{df} should be a DataFrame")
    return df


def main(provided_args: list[str] | None = None) -> None:
    """
    Main function to process command-line arguments and execute the AmpliGone tool.

    Parameters
    ----------
    provided_args : list of str, optional
        A list of command-line arguments to parse. If None, the arguments will be taken from sys.argv.

    Returns
    -------
    None

    Notes
    -----
    This function orchestrates the entire process of reading input files, processing reads, and writing output files.
    It performs the following steps:
    1. Parses the command-line arguments using the `get_args` function.
    2. Validates the provided arguments and sets the logging level.
    3. Loads the input reads and primers using concurrent futures for parallel execution.
    4. Checks the loaded reads and primers, and adjusts thread count if necessary.
    5. Processes the reads to remove primer sequences using parallel processing.
    6. Logs the results and writes the output files.

    Examples
    --------
    >>> import sys
    >>> sys.argv = ['script.py', '--input', 'input.fasta', '--primers', 'primers.fasta', '--reference', 'reference.fasta', '--output', 'output.bed']
    >>> main()
    """
    if provided_args:
        args = get_args(provided_args)
    else:
        args = get_args(sys.argv[1:])

    if len(sys.argv[1:]) < 1:
        pprint(
            f"{__prog__} was called but no arguments were given, please try again.\nUse '{__prog__} -h' to see the help document"
        )
        sys.exit(1)

    # check if verbose and quiet aren't both set
    if args.verbose is True and args.quiet is True:
        log.error(
            f"{__prog__} was given both the [green]'--verbose'[/green] and [green]'--quiet'[/green] flags. Please only use one of these flags at a time. Exiting..."
        )
        sys.exit(1)
    if args.verbose is True:
        log.setLevel("DEBUG")
        log.debug(f"Arguments: {args.__dict__}")
    if args.quiet is True:
        log.setLevel("WARNING")

    log.info(f"{__prog__} version: [blue]{__version__}[/blue]")

    # Exit if the amount of threads is None, this indicates that the system doesn't have multiple cpus available or multithreading is not enabled.
    # At least 2 threads are required for multiprocessing
    if args.threads is None:
        log.error(
            f"{__prog__} requires multithreading to run, but the system doesn't seem to have multiple cpus available (or multithreading is not enabled). Exiting..."
        )
        sys.exit(1)
    # Exit if amount of threads is less than 2, required for multiprocessing
    if args.threads < 2:
        log.error(
            f"{__prog__} requires a minimum of 2 threads for execution, but only {args.threads} thread was provided. Exiting..."
        )
        sys.exit(1)

    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        future_indexreads = pool.submit(SequenceReads, args.input)
        future_primersearch = pool.submit(
            find_or_read_primers,
            primerfile=args.primers,
            referencefile=args.reference,
            err_rate=args.error_rate,
        )
        future_scoring = pool.submit(
            AlignmentMatrix.get_scoring_matrix, args.alignment_scoring
        )

        matrix = future_scoring.result()
        primer_df = future_primersearch.result()
        indexed_reads = future_indexreads.result()

    indexed_reads = check_loaded_index(indexed_reads, args)
    primer_indexes = primer_df_to_primer_index(
        primer_df=primer_df, bind_virtual_primer=args.virtual_primers
    )
    args = check_thread_count(indexed_reads, args)

    preset: str = AlignmentPreset.get_alignment_preset(args, indexed_reads)

    args = correct_fragment_lookaround_size(args)

    log.info(
        f"Distributing {len(indexed_reads.tuples)} reads across {args.threads} threads for processing. Processing around [bold green]{round(len(indexed_reads.tuples)/args.threads)}[/bold green] reads per thread"
    )

    indexed_reads.frame = indexed_reads.frame.sample(frac=1).reset_index(drop=True)

    processed_reads = parallel_dispatcher(
        indexed_reads, args, primer_indexes, preset, matrix
    )

    total_nuc_preprocessing = sum(
        tuple(chain(indexed_reads.frame["Sequence"].str.len()))
    )
    removed_coordinates = tuple(chain(*processed_reads["Removed_coordinates"]))

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
        coord_lists_to_bed(filtered_primer_df, args.export_primers)

    processed_reads = processed_reads.drop(columns=["Removed_coordinates"])

    read_records = processed_reads.to_dict(orient="records")

    write_output(args.output, read_records, args.threads)
