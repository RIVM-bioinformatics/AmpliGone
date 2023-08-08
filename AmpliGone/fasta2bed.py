import argparse
from itertools import product

import pandas as pd
import regex as re
from Bio import Seq, SeqIO

from .func import log


def FindAmbigousOptions(seq):
    """It takes a sequence containing ambiguous nucleotides and returns a list of all possible unambiguous
    sequences

    Parameters
    ----------
    seq
        The sequence to be searched for.

    Returns
    -------
        A list of all possible combinations of the nucleotides in the sequence without ambiguous options.

    """
    ambigs = Seq.IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(ambigs.get, seq))))


def get_coords(seq, ref_seq, err_rate=0.1):
    """Finds all the positions in the reference sequence where the query sequence could be found, allowing
    for up to `err_rate` errors

    Parameters
    ----------
    seq
        The sequence you want to find in the reference sequence.
    ref_seq
        The reference sequence to search in.
    err_rate
        The maximum error rate you want to allow.

    Returns
    -------
        A set of tuples, each tuple is a match of the sequence in the reference sequence.

    """
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
    """Takes a list of sequences, and returns a list of coordinates for each sequence

    Returns
    -------
        A dataframe with the columns ref, start, end, name, score, strand, seq, revcomp

    """
    return pd.DataFrame(
        CoordListGen(*args, **kwargs),
        columns=["ref", "start", "end", "name", "score", "strand", "seq", "revcomp"],
    )


def CoordListGen(primerfile, referencefile, err_rate=0.1):
    """Takes a fasta file of primers and a fasta file of a reference sequence, and returns a list of
    dictionaries containing the coordinates of the primers in the reference sequence

    Parameters
    ----------
    primerfile
        The file containing the primers.
    referencefile
        The reference file that you want to search for the primers in.
    err_rate
        The maximum number of mismatches allowed between the primer and the reference.

    """
    keyl = ("LEFT", "PLUS", "POSITIVE", "FORWARD")
    keyr = ("RIGHT", "MINUS", "NEGATIVE", "REVERSE")

    primers = list(SeqIO.parse(primerfile, "fasta"))

    ref_file = list(SeqIO.parse(referencefile, "fasta"))
    ref_seq = [str(ref.seq) for ref in ref_file]
    ref_id = [ref.id for ref in ref_file]

    # The loop in a loop here is not a particularly efficient way of doing this.
    # But this is the easiest implementation for now, and it's not like this is a particularly
    # cpu or time intensive process anyway.
    # Might come back to this when there's more time to create a better solution.
    for ref_seq, ref_id in zip(ref_seq, ref_id):
        log.info(f"Searching for primers in reference-id: [yellow]{ref_id}[/yellow]")
        for primer in primers:
            seq = str(primer.seq)
            revcomp = Seq.reverse_complement(seq)

            coords = get_coords(seq, ref_seq, err_rate)
            rev_coords = get_coords(revcomp, ref_seq, err_rate)
            if coords and rev_coords:
                log.warning(
                    f"\tPrimer [yellow underline]{primer.id}[/yellow underline] found on both [underline]forward[/underline] and [underline]reverse strand[/underline] of [yellow underline]{ref_id}[/yellow underline].\n\t[yellow bold]Check to see if this is intended.[/yellow bold]"
                )
            if not coords and not rev_coords:
                log.warning(
                    f"\tSkipping [yellow underline]{primer.id}[/yellow underline] as it is found on neither [underline]forward[/underline] or [underline]reverse strand[/underline] of [yellow underline]{ref_id}[/yellow underline].\n\t[yellow bold]Check to see if this is intended.[/yellow bold]"
                )
                continue
            if coords and len(coords) > 1:
                log.warning(
                    f"\tPrimer [yellow underline]{primer.id}[/yellow underline] found on multiple locations on [underline]forward strand[/underline] of [yellow underline]{ref_id}[/yellow underline]: \n\t[yellow]{coords}\n\t[bold]Check to see if this is intended.[/yellow][/bold]"
                )
            if rev_coords and len(rev_coords) > 1:
                log.warning(
                    f"\tPrimer [yellow underline]{primer.id}[/yellow underline] found on multiple locations on [underline]reverse strand[/underline] of [yellow underline]{ref_id}[/yellow underline]: \n\t[yellow]{rev_coords}\n\t[bold]Check to see if this is intended.[/yellow][/bold]"
                )
            for start, end in coords.union(rev_coords):
                if any(o in primer.id for o in keyl):
                    strand = "+"
                elif any(o in primer.id for o in keyr):
                    strand = "-"
                else:
                    log.error(
                        f"Primer name {primer.id} does not contain orientation (e.g. {primer.id}_RIGHT). Consider suffixing with {keyl + keyr}"
                    )
                    continue
                yield dict(
                    ref=ref_id,
                    start=start,
                    end=end,
                    name=primer.id,
                    score=".",
                    strand=strand,
                    seq=seq,
                    revcomp=revcomp,
                )


def CoordinateListsToBed(df, outfile):
    """It takes a dataframe with columns named "ref", "start", "end", "name", "score", and "strand" and
    writes them to a file in BED (headerless tsv) format

    Parameters
    ----------
    df
        the dataframe containing the coordinates
    outfile
        the name of the file to write to

    Returns
    -------
        A dataframe with the columns ref, start, end, name, score, and strand.

    """
    return df[["ref", "start", "end", "name", "score", "strand"]].to_csv(
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
        type=float,
        help="The fraction of mismatches a primer can have with respect to the reference.",
        default=0.1,
    )

    flags = args.parse_args()

    df = MakeCoordinateLists(flags.primers, flags.reference, flags.primer_mismatch_rate)

    CoordinateListsToBed(df, flags.output)
