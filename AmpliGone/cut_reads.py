from collections import defaultdict

import mappy as mp
import pandas as pd

from .cutlery import PositionInOrAfterPrimer, PositionInOrBeforePrimer


def cut_read(
    seq,
    qual,
    PositionNeedsCutting,
    primer_list,
    position_on_reference,
    cut_direction,
    read_direction,
    cigar,
    query_start,
    query_end,
    fragment_lookaround_size,
):
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
    data,
    primer_df,
    reference,
    preset,
    scoring,
    fragment_lookaround_size,
    amplicon_type,
    workers,
):
    Frame, _threadnumber = data

    RVDict = defaultdict(set)
    FWDict = defaultdict(set)

    reference_ids = set(primer_df["ref"].unique())
    for refid in reference_ids:
        RVDict[refid] = set()
        FWDict[refid] = set()

    for _, refid, start, end, strand in primer_df[
        ["ref", "start", "end", "strand"]
    ].itertuples():

        for coord in range(start + 1, end):  # +1 because reference is 1-based
            if strand == "+":
                FWDict[refid].add(coord)
            elif strand == "-":
                RVDict[refid].add(coord)

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

    for _index, name, seq, qual in Frame[
        ["Readname", "Sequence", "Qualities"]
    ].itertuples():

        removed_coords_fw = []
        removed_coords_rv = []
        max_iter = 10  # If more iterations are needed, the sequence is discarded (not recorded)
        previous_seq = "impossible"
        cutting_is_done = False

        for i in range(max_iter):
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
                FWTuple = tuple(FWDict[hit.ctg])
                RVTuple = tuple(RVDict[hit.ctg])

                if not FWTuple or not RVTuple:
                    print(FWTuple, RVTuple, hit.ctg)

                qstart = hit.q_st
                qend = hit.q_en

                if (
                    amplicon_type == "end-to-end"
                    or (amplicon_type == "end-to-mid" and hit.strand == 1)
                    or amplicon_type == "fragmented"
                ):
                    seq, qual, removed_fw, qstart, qend = cut_read(
                        seq,
                        qual,
                        PositionNeedsCutting=PositionInOrBeforePrimer,
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
                        PositionNeedsCutting=PositionInOrAfterPrimer,
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

    ProcessedReads = pd.DataFrame(
        {
            "Readname": processed_readnames,
            "Sequence": processed_sequences,
            "Qualities": processed_qualities,
            "Removed_coordinates": removed_coords_per_read,
        }
    )

    return ProcessedReads
