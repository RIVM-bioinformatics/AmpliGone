import argparse
import os
from itertools import product
from typing import Any, Dict, Generator, Hashable, List, Set, Tuple

import pandas as pd
import regex as re
from Bio import Seq, SeqIO
from Bio.Data import IUPACData

from AmpliGone.io_ops import read_bed
from AmpliGone.log import log


def find_ambiguous_options(seq: str) -> List[str]:
    """
    Find all possible unambiguous sequences from a sequence containing ambiguous nucleotides.

    Parameters
    ----------
    seq : str
        The sequence containing ambiguous nucleotides.

    Returns
    -------
    List[str]
        A list of all possible combinations of the nucleotides in the sequence without ambiguous options.

    Examples
    --------
    >>> find_ambiguous_options('ATGCR')
    ['ATGCA', 'ATGCG']
    """
    ambigs = IUPACData.ambiguous_dna_values
    return list(map("".join, product(*[ambigs.get(nuc, nuc) for nuc in seq])))


def get_coords(seq: str, ref_seq: str, err_rate: float = 0.1) -> Set[Tuple[int, int]]:
    """
    Find all the positions in the reference sequence where the query sequence could be found, allowing
    for up to `err_rate` errors.

    Parameters
    ----------
    seq : str
        The sequence you want to find in the reference sequence.
    ref_seq : str
        The reference sequence to search in.
    err_rate : float, optional
        The maximum error rate you want to allow. Default is 0.1.

    Returns
    -------
    set of tuple of int
        A set of tuples, each tuple is a match of the sequence in the reference sequence.

    Examples
    --------
    >>> get_coords('ATCG', 'CGTATCGTACG')
    {(3, 7)}
    """
    max_errors = int(len(seq) * err_rate)

    options = find_ambiguous_options(seq.upper())

    if len(options) > 1:
        log.debug(
            f"PRIMERSEARCH :: Searching for [cyan]{len(options)}[/cyan] primer options resolved from ambiguous nucleotides in [green]{seq}[/green]"
        )

    matches = set()
    for e in range(max_errors + 1):
        for option in options:
            for match in re.finditer(
                f"(?:{str(option)}){{s<={e}}}", ref_seq, re.IGNORECASE, concurrent=True
            ):
                matches.add(match.span())
    return matches


def find_or_read_primers(
    primerfile: str, referencefile: str, err_rate: float
) -> pd.DataFrame:
    """
    Find or read primers from a given file.

    This function checks if the primer file is in BED format. If it is, it reads the file directly.
    Otherwise, it generates a DataFrame from the primer file and reference file with a specified error rate.

    Parameters
    ----------
    primerfile : str
        The path to the primer file. This can be in BED format or any other format.
    referencefile : str
        The path to the reference file. This is used when the primer file is not in BED format.
    err_rate : float
        The error rate to use when generating the DataFrame from the primer and reference files.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the primer information. The columns are ["ref", "start", "end", "name", "score", "strand", "seq", "revcomp"].

    """
    if primerfile.endswith(".bed"):
        log.info("Primer coordinates are given in BED format, skipping primer search")
        return read_bed(primerfile)
    return pd.DataFrame(
        CoordListGen(
            primerfile=primerfile, referencefile=referencefile, err_rate=err_rate
        ),
        columns=["ref", "start", "end", "name", "score", "strand", "seq", "revcomp"],
    )


def CoordListGen(
    primerfile: str, referencefile: str, err_rate: float = 0.1
) -> Generator[Dict[Hashable, Any], None, None]:
    """
    Generate a list of dictionaries containing the coordinates of primers in a reference sequence.

    Parameters
    ----------
    primerfile : str
        The name of the file containing the primers in FASTA format.
    referencefile : str
        The name of the file containing the reference sequence in FASTA format.
    err_rate : float, optional
        The maximum number of mismatches allowed between the primer and the reference. Default is 0.1.

    Yields
    ------
    dict
        A dictionary containing the coordinates of a primer in the reference sequence with keys "ref", "start", "end",
        "name", "score", "strand", "seq", and "revcomp".

    Notes
    -----
    This function takes a FASTA file of primers and a FASTA file of a reference sequence, and returns a list of
    dictionaries containing the coordinates of the primers in the reference sequence. The coordinates are returned as a
    dictionary with keys "ref", "start", "end", "name", "score", "strand", "seq", and "revcomp". The "ref" key contains
    the reference sequence ID, the "start" and "end" keys contain the start and end positions of the primer in the
    reference sequence, the "name" key contains the name of the primer, the "score" key contains a score ("."
    by default), the "strand" key contains the strand of the primer ("+" or "-"), the "seq" key contains the sequence
    of the primer, and the "revcomp" key contains the reverse complement of the primer sequence.

    Examples
    --------
    >>> for coord in CoordListGen('primers.fasta', 'reference.fasta', 0.1):
    ...     print(coord)

    """
    log.debug(f"Starting PRIMERSEARCH process\t@ ProcessID {os.getpid()}")
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
            revcomp = str(Seq.reverse_complement(seq))

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
                log.debug(
                    f"PRIMERSEARCH :: Found primer {primer.id} at coordinates {start}-{end} on {ref_id}"
                )
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


def CoordinateListsToBed(df: pd.DataFrame, outfile: str) -> None:
    """
    Write the coordinates in BED format to a file.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame containing the coordinates with columns named "ref", "start", "end", "name", "score", and "strand".
    outfile : str
        The name of the file to write to.

    Returns
    -------
    None

    Notes
    -----
    BED format is a tab-separated file format used to store genomic regions as intervals. The first three columns
    represent the reference name, start position, and end position of the interval, respectively. The remaining
    columns are optional and can be used to store additional information about the interval, such as a name, score,
    and strand.

    Examples
    --------
    >>> df = pd.DataFrame({'ref': ['chr1', 'chr2'], 'start': [100, 200], 'end': [200, 300], 'name': ['region1', 'region2'], 'score': [0.5, 0.8], 'strand': ['+', '-']})
    >>> CoordinateListsToBed(df, 'regions.bed')

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

    args.add_argument(
        "--verbose",
        action="store_true",
        help="Print debug information",
    )

    flags = args.parse_args()

    if flags.verbose:
        log.setLevel("DEBUG")

    df = find_or_read_primers(
        primerfile=flags.primers,
        referencefile=flags.reference,
        err_rate=flags.primer_mismatch_rate,
    )

    CoordinateListsToBed(df, flags.output)
