import argparse
import regex as re

from itertools import product

import pandas as pd
from Bio import SeqIO
from Bio import Seq


def FindAmbigousOptions(seq):
    ambigs = Seq.IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(ambigs.get, seq))))


def get_coords(seq, ref_seq, err_rate=0.1):
    max_errors = int(len(seq) * err_rate)

    options = FindAmbigousOptions(seq.upper())

    matches = set()
    for e in range(max_errors + 1):
        for option in options:
            for match in re.finditer(
                f"(?:{str(option)}){{s<={e}}}", ref_seq, re.IGNORECASE, concurrent=True
            ):
                matches.add(match.span())
    return matches


def MakeCoordinateLists(*args, **kwargs):
    return pd.DataFrame(CoodListGen(*args, **kwargs))


def CoodListGen(primerfile, referencefile, err_rate=0.1):
    keyl = ("LEFT", "PLUS", "POSITIVE", "FORWARD")
    keyr = ("RIGHT", "MINUS", "NEGATIVE", "REVERSE")

    primers = list(SeqIO.parse(primerfile, "fasta"))

    ref_file = SeqIO.read(referencefile, "fasta")
    ref_seq = str(ref_file.seq)

    for primer in primers:
        seq = str(primer.seq)
        revcomp = Seq.reverse_complement(seq)

        coods = get_coords(seq, ref_seq, err_rate)
        rev_coods = get_coords(revcomp, ref_seq, err_rate)
        if coods and rev_coods:
            print(
                f"Primer {primer.id} found on both forward and reverse strand of {ref_file.name}.\nCheck to see if this is intended."
            )
        if not coods and not rev_coods:
            print(
                f"Skipping {primer.id} as it is found on neither forward or reverse strand of {ref_file.name}.\nCheck to see if this is intended."
            )
            continue
        if coods and len(coods) > 1:
            print(
                f"Primer {primer.id} found on multiple locations on forward strand: \n\t{coods}.\nCheck to see if this is intended."
            )
        if rev_coods and len(rev_coods) > 1:
            print(
                f"Primer {primer.id} found on multiple locations on reverse strand: {coods}.\nCheck to see if this is intended."
            )
        for start, stop in coods.union(rev_coods):
            if any(o in primer.id for o in keyl):
                strand = "+"
            elif any(o in primer.id for o in keyr):
                strand = "-"
            else:
                print(
                    f"Primer name {primer.id} does not contain orientation (e.g. {primer.id}_RIGHT). Consider suffixing with {keyl + keyr}"
                )
                exit(1)
            yield dict(
                ref=ref_file.name,
                start=start,
                stop=stop,
                name=primer.id,
                score=".",
                strand=strand,
                seq=seq,
                revcomp=revcomp,
            )


def CoordinateListsToBed(df, outfile):
    return df[["ref", "start", "stop", "name", "score", "strand"]].to_csv(
        outfile, sep="\t", na_rep=".", header=False, index=False
    )


if __name__ == "__main__":
    import argparse

    args = argparse.ArgumentParser()

    args.add_argument(
        "--primers",
        metavar="File",
        type=str,
        help="Fasta file containing primers. Fasta headers should contain information about the orientation.",
        required=True,
    )

    args.add_argument(
        "--reference",
        metavar="File",
        type=str,
        help="Fasta file with the reference",
        required=True,
    )

    args.add_argument(
        "--output",
        metavar="File",
        type=str,
        help="The output BED file with coordinates of the primers.",
        required=True,
    )
    args.add_argument(
        "--primer-mismatch-rate",
        metavar="File",
        type=str,
        help="The fraction of mismatches a primer can have with respect to the reference.",
        default=0.1,
    )

    flags = args.parse_args()

    df = MakeCoordinateLists(
        flags.primers, flags.reference, flags["primer-mismatch-rate"]
    )

    CoordinateListsToBed(df, flags.output)
