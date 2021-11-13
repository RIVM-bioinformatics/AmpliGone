import mappy as mp
import numpy as np
import pandas as pd
import itertools

from .cutlery import (
    ReadAfterPrimer,
    ReadBeforePrimer,
    ReadInOrBeforePrimer,
    ReadInOrAfterPrimer,
    slice_fw_left,
    slice_fw_right,
    slice_rv_left,
    slice_rv_right,
)


def End_to_End(data, FWList, RVList, reference, preset, workers):
    Frame, _threadnumber = data

    readnames = Frame["Readname"].tolist()
    sequences = Frame["Sequence"].tolist()
    qualities = Frame["Qualities"].tolist()

    Aln = mp.Aligner(reference, preset=preset, best_n=1)

    processed_readnames = []
    processed_sequences = []
    processed_qualities = []
    removed_coords = []

    for i in range(len(readnames)):
        name = readnames[i]
        seq = sequences[i]
        qual = qualities[i]

        rmc = []

        looplimiter = 0
        for hit in Aln.map(seq):
            if looplimiter != 0:
                continue
            looplimiter += 1

            if hit.strand == 1:
                forward = True
            if hit.strand == -1:
                forward = False

            start = hit.r_st
            end = hit.r_en

            if forward:
                forward_modifier = 1
            else:
                forward_modifier = -1

            for FindFunction, PrimerList, ref_pos_to_cut, cut_front_modifier in (
                    [ReadInOrAfterPrimer, RVList, end, -1], # It is important that the end is cut first, as otherwise the end coordinate changes
                    [ReadInOrBeforePrimer, FWList, start, 1],
                    ):

                if forward_modifier*cut_front_modifier == 1:
                    seq_pos_to_cut = 0
                else:
                    seq_pos_to_cut = len(seq)

                cigar = hit.cigar
                if cut_front_modifier == -1:
                    cigar = reversed(cigar)

                for cigar_len, cigar_type in cigar:
                    if not FindFunction(ref_pos_to_cut, PrimerList):
                        break
                    while cigar_len > 0 and FindFunction(ref_pos_to_cut, PrimerList):
                        cigar_len -= 1

                        if cigar_type in (0,1,4):
                            seq_pos_to_cut += forward_modifier*cut_front_modifier
                        if cigar_type in (0,2,3):
                            ref_pos_to_cut += cut_front_modifier

                if forward_modifier*cut_front_modifier == 1:
                    seq = seq[seq_pos_to_cut:]
                    qual = qual[seq_pos_to_cut:]
                else:
                    seq = seq[:seq_pos_to_cut]
                    qual = qual[:seq_pos_to_cut]

        if len(seq) < 5:
            seq = np.nan
        if len(qual) < 5:
            qual = np.nan

        processed_readnames.append(name)
        processed_sequences.append(seq)
        processed_qualities.append(qual)
        removed_coords.append(rmc)

    ProcessedReads = pd.DataFrame(
        {
            "Readname": processed_readnames,
            "Sequence": processed_sequences,
            "Qualities": processed_qualities,
            "Removed_coordinates": removed_coords,
        }
    )

    ProcessedReads.dropna(subset=["Sequence", "Qualities"], inplace=True)

    return ProcessedReads


def End_to_Mid(data, FWList, RVList, reference, preset, workers):
    Frame, threadnumber = data

    readnames = Frame["Readname"].tolist()
    sequences = Frame["Sequence"].tolist()
    qualities = Frame["Qualities"].tolist()

    Aln = mp.Aligner(reference, preset=preset, best_n=1)

    processed_readnames = []
    processed_sequences = []
    processed_qualities = []
    removed_coords = []

    for i in range(len(readnames)):

        name = readnames[i]
        seq = sequences[i]
        qual = qualities[i]

        rmc = []

        looplimiter = 0
        for hit in Aln.map(seq):
            if looplimiter != 0:
                continue
            looplimiter += 1

            if hit.strand == 1:
                reverse = False
            if hit.strand == -1:
                reverse = True

            start = hit.r_st
            end = hit.r_en

            if reverse is False:

                while ReadBeforePrimer(start, FWList) is True:
                    seq, qual, start = slice_fw_left(start, seq, qual)
                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        start = hit2.r_st

                while (start in FWList) is True:
                    rmc.append(start)
                    seq, qual, start = slice_fw_left(start, seq, qual)
                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        start = hit2.r_st

            if reverse is True:

                while ReadAfterPrimer(end, RVList) is True:
                    seq, qual, end = slice_rv_right(end, seq, qual)
                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        end = hit2.r_en

                while (end in RVList) is True:
                    rmc.append(end)
                    seq, qual, end = slice_rv_right(end, seq, qual)
                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        end = hit2.r_en

            if len(seq) < 5:
                seq = np.nan
            if len(qual) < 5:
                qual = np.nan

            processed_readnames.append(name)
            processed_sequences.append(seq)
            processed_qualities.append(qual)
            removed_coords.append(rmc)

    ProcessedReads = pd.DataFrame(
        {
            "Readname": processed_readnames,
            "Sequence": processed_sequences,
            "Qualities": processed_qualities,
            "Removed_coordinates": removed_coords,
        }
    )

    ProcessedReads.dropna(subset=["Sequence", "Qualities"], inplace=True)

    return ProcessedReads
