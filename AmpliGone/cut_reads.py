import os
from collections import defaultdict
from typing import Callable, List, Tuple

import mappy as mp
import pandas as pd

from AmpliGone.log import log

from .cutlery import position_in_or_after_primer, position_in_or_before_primer


def cut_read(
    seq: str,
    qual: str,
    PositionNeedsCutting: Callable[..., bool],
    primer_list: Tuple[int, ...],
    position_on_reference: int,
    cut_direction: int,
    read_direction: int,
    cigar: List[List[int]],
    query_start: int,
    query_end: int,
    fragment_lookaround_size: int,
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
    if read_direction == cut_direction:
        # Start at the position that first matches the reference (skip soft clipped regions)
        position_on_sequence = query_start
    else:
        # End at the position that last matches the reference (skip soft clipped regions)
        position_on_sequence = query_end

    for cigar_len, cigar_type in cigar:
        while cigar_len > 0 and (
            PositionNeedsCutting(
                position_on_reference, primer_list, fragment_lookaround_size
            )
            or cigar_type not in (0, 7)  # always end with a match
        ):
            cigar_len -= 1
            removed_coords.append(position_on_reference)

            # Increment position on sequence if match/insert (in seq)/match(seq)/mismatch(seq)
            if cigar_type in (0, 1, 7, 8):
                position_on_sequence += read_direction * cut_direction

            # Increment position on reference if match/deletion (in seq)/match(seq)/mismatch(seq)
            if cigar_type in (0, 2, 7, 8):
                position_on_reference += cut_direction
        if not PositionNeedsCutting(
            position_on_reference, primer_list, fragment_lookaround_size
        ) and cigar_type in (0, 7):
            break

    if read_direction == cut_direction:
        seq = seq[position_on_sequence:]
        qual = qual[position_on_sequence:]
        query_end -= position_on_sequence
        return seq, qual, removed_coords, query_start, query_end
    seq = seq[:position_on_sequence]
    qual = qual[:position_on_sequence]
    return seq, qual, removed_coords, query_start, query_end


def CutReads(
    data: Tuple[pd.DataFrame, int],
    primer_sets: Tuple[defaultdict, defaultdict],
    reference: str,
    preset: str,
    scoring: List[int],
    fragment_lookaround_size: int,
    amplicon_type: str,
    workers: int,
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
    Frame, _threadnumber = data
    log.debug(
        f"Initiated thread {_threadnumber} @ process ID {os.getpid()} :: Processing {len(Frame)} reads."
    )

    FWDict, RVDict = primer_sets

    Aln = mp.Aligner(
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
    total_reads = len(Frame)
    for index, (_, name, seq, qual) in enumerate(
        Frame[["Readname", "Sequence", "Qualities"]].itertuples(), 1
    ):
        name: str
        seq: str
        qual: str
        if total_reads >= 10 and index % (total_reads // 10) == 0 and log.level == 10:
            completion_percentage = round(index / total_reads * 100)
            maxsize = position_in_or_before_primer.cache_info().maxsize
            currsize = position_in_or_before_primer.cache_info().currsize
            cache_usage_before = (
                currsize / maxsize * 100
                if maxsize is not None and currsize is not None
                else 0
            )
            maxsize = position_in_or_after_primer.cache_info().maxsize
            currsize = position_in_or_after_primer.cache_info().currsize
            cache_usage_after = (
                currsize / maxsize * 100
                if maxsize is not None and currsize is not None
                else 0
            )
            # todo: clean up this section of safely dividing by zero
            cache_misses = position_in_or_before_primer.cache_info().misses
            cache_hit_ratio_before = (
                (position_in_or_before_primer.cache_info().hits / cache_misses)
                if cache_misses != 0
                else 0
            )
            cache_misses = position_in_or_after_primer.cache_info().misses
            cache_hit_ratio_after = (
                (position_in_or_after_primer.cache_info().hits / cache_misses)
                if cache_misses != 0
                else 0
            )
            log.debug(
                f"Thread {_threadnumber} @ processID {os.getpid()}\t::\tReads processing {completion_percentage}% complete.\n\tMODULE {position_in_or_before_primer.__module__}.{position_in_or_before_primer.__qualname__} CACHE INFORMATION\n\t\tCache size usage = {cache_usage_before:.2f}%\n\t\tCache hit ratio = {cache_hit_ratio_before:.2f}%\n\tMODULE {position_in_or_after_primer.__module__}.{position_in_or_after_primer.__qualname__} CACHE INFORMATION\n\t\tCache size usage = {cache_usage_after:.2f}%\n\t\tCache hit ratio = {cache_hit_ratio_after:.2f}%"
            )

        removed_coords_fw = []
        removed_coords_rv = []
        previous_seq: str = "impossible"
        cutting_is_done = False

        for _ in range(max_iter):
            if cutting_is_done:
                break

            for hit in Aln.map(
                seq
            ):  # Yields only one (or no) hit, as the aligner object was initiated with best_n=1
                if len(seq) < 5 and len(qual) < 5:
                    cutting_is_done = True
                    break

                if seq == previous_seq:
                    processed_readnames.append(name)
                    processed_sequences.append(seq)
                    processed_qualities.append(qual)
                    removed_coords_per_read.append(
                        removed_coords_fw + removed_coords_rv
                    )
                    cutting_is_done = True
                    break

                previous_seq = seq

                # Fetch the primer coordinates that correspond to the reference that the read maps to
                # we're using tuples here because they are hashable

                FWTuple: Tuple[int, ...] = tuple(FWDict[hit.ctg])
                RVTuple: Tuple[int, ...] = tuple(RVDict[hit.ctg])

                if not FWTuple or not RVTuple:
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
                    seq, qual, removed_fw, qstart, qend = cut_read(
                        seq,
                        qual,
                        PositionNeedsCutting=position_in_or_before_primer,
                        primer_list=FWTuple,
                        position_on_reference=hit.r_st,
                        cut_direction=1,
                        read_direction=hit.strand,
                        cigar=hit.cigar,
                        query_start=qstart,
                        query_end=qend,
                        fragment_lookaround_size=fragment_lookaround_size,
                    )
                    removed_coords_fw.extend(removed_fw)

                if (
                    amplicon_type == "end-to-end"
                    or (amplicon_type == "end-to-mid" and hit.strand == -1)
                    or amplicon_type == "fragmented"
                ):
                    seq, qual, removed_rv, qstart, qend = cut_read(
                        seq,
                        qual,
                        PositionNeedsCutting=position_in_or_after_primer,
                        primer_list=RVTuple,
                        position_on_reference=hit.r_en,
                        cut_direction=-1,
                        read_direction=hit.strand,
                        cigar=list(reversed(hit.cigar)),
                        query_start=qstart,
                        query_end=qend,
                        fragment_lookaround_size=fragment_lookaround_size,
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
