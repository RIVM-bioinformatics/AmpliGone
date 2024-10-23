import argparse
import os
import re
from itertools import product
from typing import Dict, Generator, List, Tuple, Union

import pandas as pd
import parasail as ps
from Bio import Seq, SeqIO
from Bio.Data import IUPACData
from parasail import Cigar

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


def parse_cigar_obj(cig_obj: Cigar) -> Tuple[str, str]:
    """
    Parse a Cigar object and return the original cigar string and a cleaned cigar string.

    Parameters
    ----------
    cig_obj : Cigar
        The parasail Cigar object to parse.

    Returns
    -------
    Tuple[str, str]
        A tuple containing the original cigar string and the cleaned cigar string.

    Examples
    --------
    >>> cig_obj = Cigar("10M2D5M")
    >>> parse_cigar_obj(cig_obj)
    ('10M2D5M', '10M2D5M')

    >>> cig_obj = Cigar("2D10M2D")
    >>> parse_cigar_obj(cig_obj)
    ('2D10M2D', '10M')

    Notes
    -----
    The `parse_cigar_obj` function takes a Cigar object and extracts the original cigar string and a cleaned cigar string.
    The original cigar string represents the full cigar string as returned by the `decode` method of the Cigar object.
    The cleaned cigar string is obtained by removing any leading or trailing deletions from the original cigar string.
    This is necessary because the function uses a semi-global alignment within the parasail library.

    The function first checks if the `cigar_representation` returned by `cig_obj.decode` is a bytes object.
    If it is, the bytes object is converted to a string using the 'utf-8' encoding.
    Then, the function uses regular expressions to remove any leading or trailing deletions from the cigar string.
    The resulting original cigar string and cleaned cigar string are returned as a tuple.

    """
    cigar_representation = cig_obj.decode
    cigar_str = (
        str(cigar_representation, "utf-8")
        if isinstance(cigar_representation, bytes)
        else str(cigar_representation)
    )
    cleaned_cigar_str = re.sub(r"(^\d+D)|(\d+D$)", "", cigar_str)
    return cigar_str, cleaned_cigar_str


def count_cigar_errors(cigar: str) -> int:
    """
    Count the number of errors (insertions, deletions, and mismatches) in a CIGAR string.

    Parameters
    ----------
    cigar : str
        The CIGAR string representing the alignment.

    Returns
    -------
    int
        The total number of errors in the CIGAR string.

    Examples
    --------
    >>> count_cigar_errors("10=2I1=5D3=")
    7
    >>> count_cigar_errors("6=3X2D2=")
    5

    Notes
    -----
    This function takes a CIGAR string and counts the number of errors, which include insertions, deletions, and mismatches.
    The CIGAR string represents the alignment between two sequences, where each character in the string represents an operation.
    The "=" character represents a match, the "X" character represents a mismatch, the "I" character represents an insertion,
    and the "D" character represents a deletion.

    The function uses regular expressions to find all occurrences of insertions, deletions, and mismatches in the CIGAR string.
    It then sums the lengths of these occurrences to calculate the total number of errors.

    This function assumes that the CIGAR string is in a valid format and does not perform any error checking or validation.
    """
    insertions = sum(map(int, re.findall(r"(\d+)I", cigar)))
    deletions = sum(map(int, re.findall(r"(\d+)D", cigar)))
    mismatches = sum(map(int, re.findall(r"(\d+)X", cigar)))
    return mismatches + deletions + insertions


def get_coords(
    seq: str, ref_seq: str, err_rate: float = 0.1
) -> Tuple[str, int, int, int]:
    """
    Get the coordinates of the best primer option for a given sequence.

    Parameters
    ----------
    seq : str
        The input sequence.
    ref_seq : str
        The reference sequence.
    err_rate : float, optional
        The error rate, which determines the maximum number of errors allowed in the primer sequence, by default 0.1.

    Returns
    -------
    Tuple[str, int, int, int]
        A tuple containing the best primer option, start position, end position, and score.

    Examples
    --------
    >>> seq = "ATCG"
    >>> ref_seq = "ATCGATCG"
    >>> get_coords(seq, ref_seq, err_rate=0.1)
    ('ATCG', 0, 4, 100)

    Notes
    -----
    This function calculates the best primer option for a given sequence by allowing a maximum number of errors
    based on the error rate. It searches for primer options resolved from ambiguous nucleotides in the sequence.
    The function uses the `ps.sg_trace_scan_sat` function to perform sequence alignment and calculates the score
    and errors based on the alignment results. The start and end positions are determined from the alignment.
    The function returns the primer option with the highest score.
    """
    max_errors = int(len(seq) * err_rate)

    options = find_ambiguous_options(seq.upper())

    if len(options) > 1:
        log.debug(
            f"PRIMERSEARCH :: Searching for [cyan]{len(options)}[/cyan] primer options resolved from ambiguous nucleotides in [green]{seq}[/green]"
        )
    not_aligned_section_pattern = r"^(\d+)D"
    localresults = []
    for option in options:
        alignment = ps.sg_trace_scan_sat(
            option, ref_seq, 8, 30, ps.nuc44
        )  # TODO: validate if these gap_open and gap_extend values are appropriate for this use case. See below
        # gap_open = 8, gap_extend = 30
        # These values are chosen with the following idea: We want to allow for a few (one or two) deletions in the primer sequence vs the reference to deal with the possibility of a primer that is not perfectly matching the reference.
        # However, we only want to allow deletions of 1 nucleotide length, and we want to heavily penalize longer gaps.
        # The general idea is that a primer mismatch towards reference is okay if it compensates with one or two small deletions. However the primer should not be too different from the reference as this would indicate a bad primer design.

        cigar_str, cleaned_cigar_str = parse_cigar_obj(alignment.cigar)
        errors = count_cigar_errors(cleaned_cigar_str)
        score: int = alignment.score
        if errors > max_errors:
            score = -1

        new_start = 0
        if match := re.search(not_aligned_section_pattern, cigar_str):
            new_start = int(match[1])
        start = new_start
        end: int = alignment.end_ref + 1

        localresults.append((option, start, end, score))

    return max(localresults, key=lambda x: x[3])


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


def choose_best_fitting_coordinates(
    fw_coords: Tuple[str, int, int, int], rv_coords: Tuple[str, int, int, int]
) -> Tuple[str, int, int, int] | None:
    """
    Compares the forward and reverse coordinates and returns the best fitting coordinates based on their scores.

    Parameters
    ----------
    fw_coords : Tuple[str, int, int, int]
        The forward coordinates as a tuple of primer-sequence, start position, end position, and score.
    rv_coords : Tuple[str, int, int, int]
        The reverse coordinates as a tuple of primer-sequence, start position, end position, and score.

    Returns
    -------
    Tuple[str, int, int, int] | None
        The best fitting coordinates as a tuple of primer-sequence, start position, end position, and score, or None if no well-fitting coordinates are found.

    Examples
    --------
    >>> fw_coords = ('chr1', 100, 200, 5)
    >>> rv_coords = ('chr1', 150, 250, 7)
    >>> choose_best_fitting_coordinates(fw_coords, rv_coords)
    ('chr1', 150, 250, 7)

    >>> fw_coords = ('chr1', 100, 200, 5)
    >>> rv_coords = ('chr1', 150, 250, 3)
    >>> choose_best_fitting_coordinates(fw_coords, rv_coords)
    ('chr1', 100, 200, 5)

    Notes
    -----
    This function compares the scores of the forward and reverse coordinates.
    If the score of the forward coordinates is higher than the score of the reverse coordinates, the forward coordinates are returned.
    If the score of the reverse coordinates is higher, the reverse coordinates are returned.
    If the score of either set of coordinates is -1, that set is ignored.
    If no coordinates are found, None is returned.
    """

    # compare 'coords' and 'rev_coords'. Keep the one with the highest score, and empty the other. If both have the same score, keep both.
    if fw_coords[3] > rv_coords[3]:
        forward_coords = fw_coords
        reverse_coords = None
    elif rv_coords[3] > fw_coords[3]:
        forward_coords = None
        reverse_coords = rv_coords
    else:
        forward_coords = fw_coords
        reverse_coords = rv_coords

    # remove the record if the score is -1
    if forward_coords is not None and forward_coords[3] == -1:
        forward_coords = None
    if reverse_coords is not None and reverse_coords[3] == -1:
        reverse_coords = None

    best_fitting = forward_coords or reverse_coords
    return best_fitting or None


def CoordListGen(
    primerfile: str, referencefile: str, err_rate: float = 0.1
) -> Generator[Dict[str, Union[str, int]], None, None]:
    """
    Generate a list of coordinates for primers found in a reference sequence.

    Parameters
    ----------
    primerfile : str
        The path to the primer file in FASTA format.
    referencefile : str
        The path to the reference file in FASTA format.
    err_rate : float, optional
        The allowed error rate for primer matching. Defaults to 0.1.

    Yields
    ------
    dict
        A dictionary containing the following information for each primer:
        - ref : str
            The reference ID.
        - start : int
            The start coordinate of the primer in the reference sequence.
        - end : int
            The end coordinate of the primer in the reference sequence.
        - name : str
            The name of the primer.
        - score : int
            The alignment score of the primer.
        - strand : str
            The orientation of the primer strand ('+' or '-').
        - seq : str
            The sequence of the primer.
        - revcomp : str
            The reverse complement sequence of the primer.

    Examples
    --------
    >>> primerfile = "primers.fasta"
    >>> referencefile = "reference.fasta"
    >>> err_rate = 0.1
    >>> for primer_info in CoordListGen(primerfile, referencefile, err_rate):
    ...     print(primer_info)
    {'ref': 'reference1', 'start': 10, 'end': 20, 'name': 'primer1', 'score': 100, 'strand': '+', 'seq': 'ATCG', 'revcomp': 'CGAT'}
    {'ref': 'reference1', 'start': 30, 'end': 40, 'name': 'primer2', 'score': 90, 'strand': '-', 'seq': 'GCTA', 'revcomp': 'TAGC'}

    Notes
    -----
    This function generates a list of coordinates for primers found in a reference sequence. It takes a primer file in FASTA format and a reference file in FASTA format as input. The allowed error rate for primer matching can be specified using the `err_rate` parameter, which defaults to 0.1.
    The function uses the `get_coords` function to calculate the coordinates of the best primer option for each primer sequence. The `get_coords` function allows a maximum number of errors based on the error rate and searches for primer options resolved from ambiguous nucleotides in the sequence.
    The function yields a dictionary for each primer, containing the reference ID, start and end coordinates, name, alignment score, strand orientation, primer sequence, and reverse complement sequence.
    The examples demonstrate how to use the `CoordListGen` function to generate a list of primer coordinates for a given primer file and reference file. The yielded dictionary for each primer contains the relevant information for further analysis or processing.
    """
    log.debug(f"Starting PRIMERSEARCH process\t@ ProcessID {os.getpid()}")
    keyl = ("LEFT", "PLUS", "POSITIVE", "FORWARD")
    keyr = ("RIGHT", "MINUS", "NEGATIVE", "REVERSE")

    primers = list(SeqIO.parse(primerfile, "fasta"))

    ref_files = list(SeqIO.parse(referencefile, "fasta"))
    ref_seqs = [str(ref.seq) for ref in ref_files]
    ref_ids = [ref.id for ref in ref_files]

    # The loop in a loop here is not a particularly efficient way of doing this.
    # But this is the easiest implementation for now, and it's not like this is a particularly
    # cpu or time intensive process anyway.
    # Might come back to this when there's more time to create a better solution.
    for ref_seq, ref_id in zip(ref_seqs, ref_ids):
        log.info(f"Searching for primers in reference-id: [yellow]{ref_id}[/yellow]")
        for primer in primers:
            seq = str(primer.seq)
            revcomp = str(Seq.reverse_complement(seq))

            coords = get_coords(seq, ref_seq, err_rate)
            rev_coords = get_coords(revcomp, ref_seq, err_rate)

            best_fitting = choose_best_fitting_coordinates(coords, rev_coords)

            if not best_fitting:
                log.warning(
                    f"\tSkipping [yellow underline]{primer.id}[/yellow underline] as it is found on neither [underline]forward[/underline] or [underline]reverse strand[/underline] of [yellow underline]{ref_id}[/yellow underline].\n\t[yellow bold]Check to see if this is intended.[/yellow bold]"
                )
                continue

            _, start, end, score = best_fitting

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
                f"PRIMERSEARCH :: Found primer [yellow]{primer.id}[/yellow] at coordinates [cyan]{start}-{end}[/cyan] with alignment-score [cyan]{score}[/cyan] on [yellow]{ref_id}[/yellow]"
            )
            yield {
                "ref": ref_id,
                "start": start,
                "end": end,
                "name": primer.id,
                "score": score,
                "strand": strand,
                "seq": seq,
                "revcomp": revcomp,
            }


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


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--primers",
        metavar="File",
        type=str,
        help="Fasta file containing primers. Fasta headers should contain information about the orientation.",
        required=True,
    )

    parser.add_argument(
        "--reference",
        metavar="File",
        type=str,
        help="Fasta file with the reference",
        required=True,
    )

    parser.add_argument(
        "--output",
        metavar="File",
        type=str,
        help="The output BED file with coordinates of the primers.",
        required=True,
    )
    parser.add_argument(
        "--primer-mismatch-rate",
        metavar="File",
        type=float,
        help="The fraction of mismatches a primer can have with respect to the reference.",
        default=0.1,
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print debug information",
    )

    return parser.parse_args(args)


def main(args: list[str] | None = None) -> None:
    flags = parse_args(args)
    if flags.verbose:
        log.setLevel("DEBUG")

    df = find_or_read_primers(
        primerfile=flags.primers,
        referencefile=flags.reference,
        err_rate=flags.primer_mismatch_rate,
    )

    CoordinateListsToBed(df, flags.output)


if __name__ == "__main__":
    main()
