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
    Frame, threadnumber = data
    FWSet, RVSet = set(FWList), set(RVList)

    readnames = Frame["Readname"].tolist()
    sequences = Frame["Sequence"].tolist()
    qualities = Frame["Qualities"].tolist()

    Aln = mp.Aligner(reference, preset=preset, best_n=1)

    processed_readnames = []
    processed_sequences = []
    processed_qualities = []
    removed_coords = []

    new_method = 0

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
                if ReadInOrAfterPrimer(end, RVList):
                    seq_pos_to_cut = len(seq)
                    ref_pos_to_cut = end
                    for hit2 in Aln.map(seq):
                        for cigar_len, cigar_type in reversed(hit2.cigar):
                            if not ReadInOrAfterPrimer(ref_pos_to_cut, RVList):
                                break
                            while cigar_len > 0 and ReadInOrAfterPrimer(ref_pos_to_cut, RVList):
                                cigar_len -= 1
                                if cigar_type == 0 or cigar_type == 3: # MATCH or MISMATCH
                                    seq_pos_to_cut -= 1
                                    ref_pos_to_cut -= 1
                                elif cigar_type == 1: # Insertion in seq
                                    seq_pos_to_cut -= 1
                                elif cigar_type == 2: # Deletion in seq
                                    ref_pos_to_cut -= 1
                                else:
                                    raise ValueError(f'Error in finding position to cut read:\nCigar type {cigar_type} does not exist')
                        seq = seq[:seq_pos_to_cut]
                        qual = qual[:seq_pos_to_cut]
                        new_method += 1

                if ReadInOrBeforePrimer(start, FWList):
                    seq_pos_to_cut = 0
                    ref_pos_to_cut = start
                    for cigar_len, cigar_type in hit.cigar:
                        if not ReadInOrBeforePrimer(ref_pos_to_cut, FWList):
                            break
                        while cigar_len > 0 and ReadInOrBeforePrimer(ref_pos_to_cut, FWList):
                            cigar_len -= 1
                            if cigar_type == 0 or cigar_type == 3: # MATCH or MISMATCH
                                seq_pos_to_cut += 1
                                ref_pos_to_cut += 1
                            elif cigar_type == 1: # Insertion in seq
                                seq_pos_to_cut += 1
                            elif cigar_type == 2: # Deletion in seq
                                ref_pos_to_cut += 1
                            else:
                                raise ValueError(f'Error in finding position to cut read:\nCigar type {cigar_type} does not exist')
                    seq = seq[seq_pos_to_cut:]
                    qual = qual[seq_pos_to_cut:]
                    new_method += 1

            if reverse is True:

                if ReadInOrBeforePrimer(start, FWList):
                    seq_pos_to_cut = len(seq)
                    ref_pos_to_cut = start
                    for hit2 in Aln.map(seq):
                        for cigar_len, cigar_type in hit2.cigar:
                            if not ReadInOrBeforePrimer(ref_pos_to_cut, FWList):
                                break
                            while cigar_len > 0 and ReadInOrBeforePrimer(ref_pos_to_cut, FWList):
                                cigar_len -= 1
                                if cigar_type == 0 or cigar_type == 3: # MATCH or MISMATCH
                                    seq_pos_to_cut -= 1
                                    ref_pos_to_cut += 1
                                elif cigar_type == 1: # Insertion in seq
                                    seq_pos_to_cut -= 1
                                elif cigar_type == 2: # Deletion in seq
                                    ref_pos_to_cut += 1
                                else:
                                    raise ValueError(f'Error in finding position to cut read:\nCigar type {cigar_type} does not exist')
                        seq = seq[:seq_pos_to_cut]
                        qual = qual[:seq_pos_to_cut]
                        new_method += 1

                if ReadInOrAfterPrimer(end, RVList):
                    seq_pos_to_cut = 0
                    ref_pos_to_cut = end
                    for hit2 in Aln.map(seq):
                        for cigar_len, cigar_type in hit2.cigar:
                            if not ReadInOrAfterPrimer(ref_pos_to_cut, RVList):
                                break
                            while cigar_len > 0 and ReadInOrAfterPrimer(ref_pos_to_cut, RVList):
                                cigar_len -= 1
                                if cigar_type == 0 or cigar_type == 3: # MATCH or MISMATCH
                                    seq_pos_to_cut += 1
                                    ref_pos_to_cut -= 1
                                elif cigar_type == 1: # Insertion in seq
                                    seq_pos_to_cut += 1
                                elif cigar_type == 2: # Deletion in seq
                                    ref_pos_to_cut -= 1
                                else:
                                    raise ValueError(f'Error in finding position to cut read:\nCigar type {cigar_type} does not exist')
                        seq = seq[seq_pos_to_cut:]
                        qual = qual[seq_pos_to_cut:]
                        new_method += 1

                # while ReadAfterPrimer(end, RVList) is True:
                #     seq, qual, end = slice_rv_right(end, seq, qual)
                #     hitlimit = 0
                #     for hit2 in Aln.map(seq):
                #         if hitlimit != 0:
                #             continue
                #         hitlimit += 1
                #         end = hit2.r_en

                # while (end in RVSet) is True:
                #     rmc.append(end)
                #     seq, qual, end = slice_rv_right(end, seq, qual)
                #     hitlimit = 0
                #     for hit2 in Aln.map(seq):
                #         if hitlimit != 0:
                #             continue
                #         hitlimit += 1
                #         end = hit2.r_en

        if len(seq) < 5:
            seq = np.nan
        if len(qual) < 5:
            qual = np.nan

        processed_readnames.append(name)
        processed_sequences.append(seq)
        processed_qualities.append(qual)
        removed_coords.append(rmc)

    print(f"Removed with new_method: {new_method}")
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
