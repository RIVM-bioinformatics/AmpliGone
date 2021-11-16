import mappy as mp
import numpy as np
import pandas as pd

from .cutlery import (
    ReadAfterPrimer,
    ReadBeforePrimer,
    PositionInOrBeforePrimer,
    PositionInOrAfterPrimer,
    slice_fw_left,
    slice_rv_right,
)

def cut_read(seq, qual, PositionNeedsCutting, primer_list, position_on_reference, cut_direction, read_direction, cigar):
    removed_coords = []

    # Whether to start at the end or at the start of the read sequence
    if read_direction == cut_direction:
        position_on_sequence = 0
    else:
        position_on_sequence = len(seq)

    for cigar_len, cigar_type in cigar:
        if not PositionNeedsCutting(position_on_reference, primer_list):
            break

        while cigar_len > 0 and PositionNeedsCutting(position_on_reference, primer_list):
            cigar_len -= 1
            removed_coords.append(position_on_reference)

            # Increment position on sequence if match/insert(in seq)/mismatch
            if cigar_type in (0,1,4):
                position_on_sequence += read_direction*cut_direction

            # Increment position on reference if match/insert(in seq)/mismatch
            if cigar_type in (0,2,3):
                position_on_reference += cut_direction

    if read_direction == cut_direction:
        seq = seq[position_on_sequence:]
        qual = qual[position_on_sequence:]
    else:
        seq = seq[:position_on_sequence]
        qual = qual[:position_on_sequence]
    return seq,qual,removed_coords

def CutReads(data, FWList, RVList, reference, preset, workers, amplicon_type):
    Frame, _threadnumber = data

    readnames = Frame["Readname"].tolist()
    sequences = Frame["Sequence"].tolist()
    qualities = Frame["Qualities"].tolist()

    Aln = mp.Aligner(reference, preset=preset, best_n=1)

    processed_readnames = []
    processed_sequences = []
    processed_qualities = []
    removed_coords_per_read = [] # A list of lists

    for name, seq, qual in zip(readnames, sequences, qualities):

        for hit in Aln.map(seq): # Yields only one hit, as the aligner object was initiated with best_n=1

            removed_coords_fw = removed_coords_rv = []

            if amplicon_type == 'end-to-end' or (amplicon_type == 'end_to_mid' and hit.strand == 1):
                seq, qual, removed_coords_fw = cut_read(
                    seq, qual,
                    PositionNeedsCutting=PositionInOrBeforePrimer,
                    primer_list=FWList,
                    position_on_reference=hit.r_st,
                    cut_direction=1,
                    read_direction=hit.strand,
                    cigar=hit.cigar,
                    )

            if amplicon_type == 'end-to-end' or (amplicon_type == 'end_to_mid' and hit.strand == -1):
                seq, qual, removed_coords_rv = cut_read(
                    seq, qual,
                    PositionNeedsCutting=PositionInOrAfterPrimer,
                    primer_list=RVList,
                    position_on_reference=hit.r_en,
                    cut_direction=-1,
                    read_direction=hit.strand,
                    cigar=list(reversed(hit.cigar)),
                    )


        if len(seq) >= 5 and len(qual) >= 5:
            processed_readnames.append(name)
            processed_sequences.append(seq)
            processed_qualities.append(qual)
            removed_coords_per_read.append(removed_coords_fw + removed_coords_rv)

    ProcessedReads = pd.DataFrame(
        {
            "Readname": processed_readnames,
            "Sequence": processed_sequences,
            "Qualities": processed_qualities,
            "Removed_coordinates": removed_coords_per_read,
        }
    )

    return ProcessedReads
