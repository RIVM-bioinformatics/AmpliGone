import mappy as mp
import pandas as pd
import numpy as np
from tqdm import tqdm
from .cutlery import ReadBeforePrimer, ReadAfterPrimer
from .cutlery import slice_fw_left, slice_fw_right, slice_rv_left, slice_rv_right


def End_to_End(data, FWList, RVList, reference, preset, workers):
    Frame, threadnumber = data

    readnames = Frame["Readname"].tolist()
    sequences = Frame["Sequence"].tolist()
    qualities = Frame["Qualities"].tolist()

    Aln = mp.Aligner(reference, preset=preset, best_n=1)

    processed_readnames = []
    processed_sequences = []
    processed_qualities = []

    for i in range(len(readnames)):
        name = readnames[i]
        seq = sequences[i]
        qual = qualities[i]

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

                while ReadAfterPrimer(end, RVList) is True:
                    seq, qual, end = slice_fw_right(end, seq, qual)

                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        end = hit2.r_en

                while (end in RVList) is True:
                    seq, qual, end = slice_fw_right(end, seq, qual)
                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        end = hit2.r_en

                while (start in FWList) is True:
                    seq, qual, start = slice_fw_left(start, seq, qual)
                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        start = hit2.r_st

            if reverse is True:

                while ReadBeforePrimer(start, FWList) is True:
                    seq, qual, start = slice_rv_left(start, seq, qual)
                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        start = hit2.r_st

                while ReadAfterPrimer(end, RVList) is True:
                    seq, qual, end = slice_rv_right(end, seq, qual)
                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        end = hit2.r_en

                while (start in FWList) is True:
                    seq, qual, start = slice_rv_left(start, seq, qual)
                    hitlimit = 0
                    for hit2 in Aln.map(seq):
                        if hitlimit != 0:
                            continue
                        hitlimit += 1
                        start = hit2.r_st

                while (end in RVList) is True:
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

    ProcessedReads = pd.DataFrame(
        {
            "Readname": processed_readnames,
            "Sequence": processed_sequences,
            "Qualities": processed_qualities,
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

    for i in range(len(readnames)):

        name = readnames[i]
        seq = sequences[i]
        qual = qualities[i]

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

    ProcessedReads = pd.DataFrame(
        {
            "Readname": processed_readnames,
            "Sequence": processed_sequences,
            "Qualities": processed_qualities,
        }
    )

    ProcessedReads.dropna(subset=["Sequence", "Qualities"], inplace=True)

    return ProcessedReads
