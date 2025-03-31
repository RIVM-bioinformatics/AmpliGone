"""
This module provides functions to cut read sequences based on primer locations and reference mapping. It includes
functions to handle sequence reads, quality strings, and CIGAR information to accurately cut and process reads.

Functions
---------
cut_read(seq: str, qual: str, position_needs_cutting: Callable[..., bool], primer_list: Tuple[int, ...], position_on_reference: int, cut_direction: int, read_direction: int, cigar: List[List[int]], query_start: int, query_end: int, fragment_lookaround_size: int) -> Tuple[str, str, List[int], int, int]
    Cut a read sequence and quality string based on read-direction, CIGAR information, orientation, and fragment position.

cut_reads(data: Tuple[pd.DataFrame, int], primer_sets: Tuple[defaultdict, defaultdict], reference: str, preset: str, scoring: List[int], fragment_lookaround_size: int, amplicon_type: str) -> pd.DataFrame
    Cut reads based on primer locations and reference mapping.

Notes
-----
- The `cut_read` function processes individual read sequences and quality strings, cutting them based on the provided
  parameters.
- The `cut_reads` function processes a DataFrame of reads, cutting them based on primer locations and reference mapping.
- The module uses the `mappy` library for sequence alignment and the `pandas` library for data manipulation.
- The `position_in_or_before_primer` and `position_in_or_after_primer` functions are used to determine if a position
  needs cutting based on primer locations.
"""

# pylint: disable=E1120
# pylint doesnt understand that the position_in_or_after_primer function is not being called the lines around 174
import os
from collections import defaultdict
from dataclasses import dataclass
from functools import _lru_cache_wrapper
from typing import Callable, List, Tuple

# mappy is a C extension, so it is added to the pylint extension allow list
import mappy as mp
import pandas as pd

from AmpliGone.log import log

from .cutlery import position_in_or_after_primer, position_in_or_before_primer


@dataclass
class Read:
    """
    A class to represent a sequencing read.

    Attributes
    ----------
    name : str
        The name or identifier of the read.
    seq : str
        The nucleotide sequence of the read.
    qual : str
        The quality scores of the read, encoded as a string.
    """

    name: str
    seq: str
    qual: str


@dataclass
class CuttingParameters:
    """
    A class to represent the parameters required for cutting reads.

    Attributes
    ----------
    position_needs_cutting : Callable[..., bool]
        A function that determines if a position needs cutting based on various criteria.

    primer_list : Tuple[int, ...]
        A tuple containing the list of primer positions.

    position_on_reference : int
        The position on the reference sequence.

    cut_direction : int
        The direction of the cut (e.g., 1 for forward, -1 for reverse).

    read_direction : int
        The direction of the read (e.g., 1 for forward, -1 for reverse).

    cigar : List[List[int]]
        The CIGAR string representing the alignment, as a list of operations.

    query_range : dict[str, int]
        A dictionary containing the start and end positions of the query sequence.

    fragment_lookaround_size : int
        The size of the fragment lookaround region.
    """

    position_needs_cutting: Callable[..., bool]
    primer_list: Tuple[int, ...]
    position_on_reference: int
    cut_direction: int
    read_direction: int
    cigar: List[List[int]]
    query_range: dict[str, int]
    fragment_lookaround_size: int


def cut_read(
    read: Read, params: CuttingParameters
) -> Tuple[str, str, List[int], int, int]:
    """
    Cut a read sequence and quality string based read-direction, cigar-information, orientation, and fragment position.

    Parameters
    ----------
    seq : str
        The read sequence.

    qual : str
        The quality string.

    PositionNeedsCutting : Callable[..., bool]
        A function that returns True if the position on the reference needs to be cut.

    primer_list : Tuple[int, ...]
        A tuple of integers representing the positions of primers on the reference.

    position_on_reference : int
        The position on the reference where the read sequence starts.

    cut_direction : int
        The direction in which the read sequence needs to be cut.

    read_direction : int
        The direction in which the read sequence is read.

    cigar : List[List[int]]
        A list of lists representing the CIGAR string.

    query_start : int
        The start position of the read sequence on the query.

    query_end : int
        The end position of the read sequence on the query.

    fragment_lookaround_size : int
        The size of the fragment lookaround.

    Returns
    -------
    Tuple[str, str, List[int], int, int]
        A tuple containing the cut read sequence, the cut quality string, a list of removed coordinates,
        the new query start position, and the new query end position.
    """
    removed_coords = []

    # Whether to start at the end or at the start of the read sequence
    if params.read_direction == params.cut_direction:
        # Start at the position that first matches the reference (skip soft clipped regions)
        position_on_sequence = params.query_range["start"]
    else:
        # End at the position that last matches the reference (skip soft clipped regions)
        position_on_sequence = params.query_range["end"]

    for cigar_len, cigar_type in params.cigar:
        while cigar_len > 0 and (
            params.position_needs_cutting(
                params.position_on_reference,
                params.primer_list,
                params.fragment_lookaround_size,
            )
            or cigar_type not in (0, 7)  # always end with a match
        ):
            cigar_len -= 1
            removed_coords.append(params.position_on_reference)

            # Increment position on sequence if match/insert (in seq)/match(seq)/mismatch(seq)
            if cigar_type in (0, 1, 7, 8):
                position_on_sequence += params.read_direction * params.cut_direction

            # Increment position on reference if match/deletion (in seq)/match(seq)/mismatch(seq)
            if cigar_type in (0, 2, 7, 8):
                params.position_on_reference += params.cut_direction
        if not params.position_needs_cutting(
            params.position_on_reference,
            params.primer_list,
            params.fragment_lookaround_size,
        ) and cigar_type in (0, 7):
            break

    if params.read_direction == params.cut_direction:
        read.seq = read.seq[position_on_sequence:]
        read.qual = read.qual[position_on_sequence:]
        params.query_range["end"] -= position_on_sequence
        return (
            read.seq,
            read.qual,
            removed_coords,
            params.query_range["start"],
            params.query_range["end"],
        )
    read.seq = read.seq[:position_on_sequence]
    read.qual = read.qual[:position_on_sequence]
    return (
        read.seq,
        read.qual,
        removed_coords,
        params.query_range["start"],
        params.query_range["end"],
    )


def log_cache_info(index: int, total_reads: int, _threadnumber: int) -> None:
    """
    Logs cache information for the primer position functions.

    Parameters
    ----------
    index : int
        The current index of the read being processed.

    total_reads : int
        The total number of reads to be processed.

    _threadnumber : int
        The thread number of the current process.

    Returns
    -------
    None
        This function does not return any value. It logs the cache information.

    Notes
    -----
    This function logs the completion percentage of read processing and the cache usage and hit ratio
    for the `position_in_or_before_primer` and `position_in_or_after_primer` functions. It also handles
    potential division by zero errors when calculating cache hit ratios.
    """

    def _get_cache_stats(func: _lru_cache_wrapper) -> Tuple[int, int, int]:
        """
        Get cache statistics for a given function.

        Parameters
        ----------
        func : Callable[..., bool]
            The function to get cache statistics for.

        Returns
        -------
        Tuple[int, int, int]
            A tuple containing the current size of the cache, the number of cache misses, and the number of cache hits.
        """
        cache_info = func.cache_info()
        currsize = cache_info.currsize
        misses = cache_info.misses
        hits = cache_info.hits
        return currsize, misses, hits

    def _get_hit_ratio(hits: int, misses: int) -> float:
        """
        Calculates the hit ratio.

        Parameters
        ----------
        hits : int
            The number of hits.
        misses : int
            The number of misses.

        Returns
        -------
        float
            The hit ratio, as a percentage.
            Returns 0 if either misses or hits are 0 to avoid division by zero.
        """
        return (hits / (misses + hits)) * 100 if misses != 0 and hits != 0 else 0

    completion_percentage = round(index / total_reads * 100)

    before_usedsize, before_misses, before_hits = _get_cache_stats(
        position_in_or_before_primer
    )
    after_usedsize, after_misses, after_hits = _get_cache_stats(
        position_in_or_after_primer
    )

    hit_ratio_before = _get_hit_ratio(before_hits, before_misses)
    hit_ratio_after = _get_hit_ratio(after_hits, after_misses)

    log.debug(
        # mypy doesnt understand that the position_in_or_before_primer has a __qualname__ attribute,
        # because it thinks its the wrapper (lru_cache) function, which does not have a __qualname__ attribute
        f"Thread {_threadnumber} @ processID {os.getpid()}\t::\t"
        f"Reads processing {completion_percentage}% complete.\n\t"
        f"MODULE {position_in_or_before_primer.__module__}.{position_in_or_before_primer.__qualname__} "  # type: ignore[attr-defined]
        f"CACHE INFORMATION\n\t\t{before_usedsize} unique records stored in cache\n\t\t"
        f"{hit_ratio_before:.2f}% cache hit ratio\n\t"
        f"MODULE {position_in_or_after_primer.__module__}.{position_in_or_after_primer.__qualname__} "
        f"CACHE INFORMATION\n\t\t{after_usedsize} unique records stored in cache\n\t\t"
        f"{hit_ratio_after:.2f}% cache hit ratio"
    )


def cut_reads(
    data: Tuple[pd.DataFrame, int],
    primer_sets: Tuple[defaultdict, defaultdict],
    reference: str,
    preset: str,
    scoring: List[int],
    fragment_lookaround_size: int,
    amplicon_type: str,
) -> pd.DataFrame:
    """
    Cut reads based on primer locations and reference mapping.

    Parameters
    ----------
    data : Tuple[pd.DataFrame, int]
        A tuple containing a pandas DataFrame with columns "Readname", "Sequence", and "Qualities",
        and an integer representing the thread number.

    primer_sets : Tuple[defaultdict, defaultdict]
        A tuple containing two defaultdicts, one for forward primers and one for reverse primers. These defaultdicts contain the primer coordinates to remove.

    reference : str
        The reference genome sequence.

    preset : str
        The preset used for minimap2 alignment.

    scoring : List[int]
        The scoring matrix used for minimap2 alignment.

    fragment_lookaround_size : int
        The number of bases to look around a fragment when cutting reads.

    amplicon_type : str
        The type of amplicon, either "end-to-end", "end-to-mid", or "fragmented".

    workers : int
        The number of workers to use for parallel processing.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame with columns "Readname", "Sequence", "Qualities", and "Removed_coordinates",
        representing the processed reads and the coordinates that were removed.
    """
    frame, _threadnumber = data
    log.debug(
        f"Initiated thread {_threadnumber} @ process ID {os.getpid()} :: Processing {len(frame)} reads."
    )

    fw_dict, rv_dict = primer_sets

    aligner = mp.Aligner(
        reference,
        preset=preset,
        best_n=1,
        scoring=scoring,
        extra_flags=0x4000000,  # Distinguish between match and mismatch: MM_F_EQX flag in minimap2
    )
    processed_readnames = []
    processed_sequences = []
    processed_qualities = []
    removed_coords_per_read = []  # A list of lists

    max_iter = (
        10  # If more iterations are needed, the sequence is discarded (not recorded)
    )
    total_reads = len(frame)
    for index, (_, name, seq, qual) in enumerate(
        frame[["Readname", "Sequence", "Qualities"]].itertuples(), 1
    ):
        if (
            total_reads >= 10 and index % (total_reads // 10) == 0 and log.level == 10
        ):  # if the log level is set to DEBUG, print the cache info every 10% of the reads.
            # This may be a bit overwhelming, but the cache info could be very informative to find potential performance issues.
            log_cache_info(index, total_reads, _threadnumber)

        if len(seq) < 42:
            # Length of the read has to be at least ~42bp because the default k-mer size for the short reads preset (sr) is 21.

            log.debug(f"Read with name '{name}' is too short to be processed.")
            continue

        read = Read(name, seq, qual)
        removed_coords_fw: list[int | None] = []
        removed_coords_rv: list[int | None] = []
        previous_seq: str = "impossible"
        cutting_is_done = False

        for _ in range(max_iter):
            if cutting_is_done:
                break

            for hit in aligner.map(
                read.seq
            ):  # Yields only one (or no) hit, as the aligner object was initiated with best_n=1
                if len(read.seq) < 5 and len(read.qual) < 5:
                    cutting_is_done = True
                    break

                if read.seq == previous_seq:
                    processed_readnames.append(read.name)
                    processed_sequences.append(read.seq)
                    processed_qualities.append(read.qual)
                    removed_coords_per_read.append(
                        removed_coords_fw + removed_coords_rv
                    )
                    cutting_is_done = True
                    break

                previous_seq = read.seq

                # Fetch the primer coordinates that correspond to the reference that the read maps to
                # we're using tuples here because they are hashable

                fw_tuple: Tuple[int, ...] = tuple(fw_dict[hit.ctg])
                rv_tuple: Tuple[int, ...] = tuple(rv_dict[hit.ctg])

                if not fw_tuple or not rv_tuple:
                    log.debug(
                        f"Thread {_threadnumber} @ processID {os.getpid()}\t::\tRead with name '{name}' aligns to '{hit.ctg}', but there are no primers affiliated with '{hit.ctg}'."
                    )
                    continue

                qstart: int = hit.q_st
                qend: int = hit.q_en

                if (
                    amplicon_type == "end-to-end"
                    or (amplicon_type == "end-to-mid" and hit.strand == 1)
                    or amplicon_type == "fragmented"
                ):
                    params = CuttingParameters(
                        position_needs_cutting=position_in_or_before_primer,
                        primer_list=fw_tuple,
                        position_on_reference=hit.r_st,
                        cut_direction=1,
                        read_direction=hit.strand,
                        cigar=hit.cigar,
                        query_range={"start": qstart, "end": qend},
                        fragment_lookaround_size=fragment_lookaround_size,
                    )
                    read.seq, read.qual, removed_fw, qstart, qend = cut_read(
                        read, params
                    )
                    removed_coords_fw.extend(removed_fw)

                if (
                    amplicon_type == "end-to-end"
                    or (amplicon_type == "end-to-mid" and hit.strand == -1)
                    or amplicon_type == "fragmented"
                ):
                    params = CuttingParameters(
                        position_needs_cutting=position_in_or_after_primer,
                        primer_list=rv_tuple,
                        position_on_reference=hit.r_en,
                        cut_direction=-1,
                        read_direction=hit.strand,
                        cigar=list(reversed(hit.cigar)),
                        query_range={"start": qstart, "end": qend},
                        fragment_lookaround_size=fragment_lookaround_size,
                    )
                    read.seq, read.qual, removed_rv, qstart, qend = cut_read(
                        read, params
                    )
                    removed_coords_rv.extend(removed_rv)

    return pd.DataFrame(
        {
            "Readname": processed_readnames,
            "Sequence": processed_sequences,
            "Qualities": processed_qualities,
            "Removed_coordinates": removed_coords_per_read,
        }
    )
