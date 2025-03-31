"""
This module provides functions and classes for determining the optimal alignment preset for sequencing reads.

Functions
---------
get_alignment_preset(input_args: argparse.Namespace, indexed_reads: SequenceReads) -> str
    Get the alignment preset for the given reads.

find_preset(threads: int, data: pd.DataFrame) -> str
    Find the preset based on the statistics calculated from the input data.

_qual_to_ord_dispatcher(qdata: str, threads: int) -> List[int]
    Convert a string of characters to a list of ASCII values minus 33 using multiple threads.

_process_chunk(chunk: str) -> List[int]
    Converts each character in the given chunk to its corresponding ASCII value minus 33.

_determine_preset(avg_len: float, avg_qual: float, quality_range: int, length_range: int) -> str
    Determine the preset based on the calculated statistics.

_calc_avg_read_length(sequence_list: List[str]) -> float
    Calculate the average length of a list of sequences.

_calc_avg_read_qual(quality_list: List[int]) -> float
    Calculate the average quality score of a list of quality scores.

_get_unique_quality_scores(quality_list: List[int]) -> int
    Return the number of unique quality scores in a list.

_get_unique_read_lengths(sequence_list: List[str]) -> int
    Return the number of unique lengths of strings in a list.

_sequence_statistics_dispatcher(reads_list: List[str], ordinal_qualities_list: List[int], threads: int) -> Tuple[float, float, int, int]
    Calculate sequence statistics using multiple threads.

Notes
-----
This module is designed to handle sequencing reads and determine the optimal alignment preset based on the quality and length of the reads. It uses parallel processing to efficiently handle large datasets and calculate necessary statistics.

Examples
--------
>>> import pandas as pd
>>> data = pd.DataFrame({'Sequence': ['ATCG', 'GCTA'], 'Qualities': ['!@#$%', '&*()']})
>>> find_preset(4, data)
'sr'
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
from typing import Generator, List, Tuple

import pandas as pd

from AmpliGone.io_ops import SequenceReads
from AmpliGone.log import log


def get_alignment_preset(
    input_args: argparse.Namespace, indexed_reads: SequenceReads
) -> str:
    """
    Get the alignment preset for the given reads.

    Parameters
    ----------
    input_args : argparse.Namespace
        The input arguments.
    indexed_reads : LoadData
        The indexed reads.

    Returns
    -------
    str
        The alignment preset.

    Raises
    ------
    None

    Notes
    -----
    If the `alignment_preset` argument is provided in `input_args`, it will be returned directly.
    Otherwise, the function will find the optimal alignment preset for the given reads.

    The optimal alignment preset is determined by sampling a subset of the reads and using the
    `find_preset` function to calculate the preset based on the sampled reads.

    The sample size is determined by taking the minimum of the total number of reads and 15000.

    """
    if input_args.alignment_preset is not None:
        # this check is mostly for mypy to understand that alignment_preset is not Any
        if not isinstance(input_args.alignment_preset, str):
            raise TypeError("alignment_preset should be a string")
        return input_args.alignment_preset
    log.info("Finding optimal alignment-preset for the given reads")
    sample_size = min(len(indexed_reads.tuples), 15000)
    return find_preset(input_args.threads, indexed_reads.frame.sample(n=sample_size))


def find_preset(threads: int, data: pd.DataFrame) -> str:
    """
    Find the preset based on the statistics calculated from the input data.

    Parameters
    ----------
    threads : int
        Number of threads to be used for processing the data.
    data : pd.DataFrame
        Pandas DataFrame containing at least two columns: "Sequence" and "Qualities".
        The "Sequence" column should contain strings representing sequences, and the "Qualities" column
        should contain strings representing quality data associated with the sequences.

    Returns
    -------
    str
        The preset determined based on the calculated statistics.

    Notes
    -----
    This function takes in a number of threads and a DataFrame, extracts read data, calculates statistics,
    and determines a preset based on the statistics.

    The preset is determined by calling the `_determine_preset` function with the calculated statistics as arguments.

    See Also
    --------
    _determine_preset : Function to determine the preset based on the calculated statistics.

    Examples
    --------
    >>> import pandas as pd
    >>> data = pd.DataFrame({'Sequence': ['ATCG', 'GCTA'], 'Qualities': ['!@#$%', '&*()']})
    >>> find_preset(4, data)
    'sr'
    """

    def _extract_read_data(data: pd.DataFrame) -> Tuple[List[str], str]:
        """
        Extract read data from the input DataFrame.

        Parameters
        ----------
        data : pd.DataFrame
            Pandas DataFrame containing at least two columns: "Sequence" and "Qualities".

        Returns
        -------
        Tuple[List[str], str]
            A tuple containing two elements:
            1. `ReadList` : List[str]
                A list of strings extracted from the "Sequence" column of the input DataFrame.
            2. `QualData` : str
                A single string obtained by joining all the elements in the "Qualities" column of the input DataFrame.

        Notes
        -----
        This function takes a pandas DataFrame as input and returns a tuple containing a list of sequences
        and a concatenated string of qualities.

        Examples
        --------
        >>> import pandas as pd
        >>> data = pd.DataFrame({'Sequence': ['ATCG', 'GCTA'], 'Qualities': ['!@#$%', '&*()']})
        >>> _extract_read_data(data)
        (['ATCG', 'GCTA'], '!@#$%&*()')
        """
        list_of_reads: List[str] = data["Sequence"].tolist()
        qualities_str: str = "".join(data["Qualities"].tolist())
        return list_of_reads, qualities_str

    reads_list, quality_data = _extract_read_data(data)
    ord_quality_list: list[int] = _qual_to_ord_dispatcher(quality_data, threads)
    avg_len, avg_qual, quality_range, length_range = _sequence_statistics_dispatcher(
        reads_list, ord_quality_list, threads
    )
    return _determine_preset(avg_len, avg_qual, quality_range, length_range)


def _qual_to_ord_dispatcher(qdata: str, threads: int) -> list[int]:
    """
    Convert a string of characters to a list of ASCII values minus 33 using multiple threads.

    Parameters
    ----------
    qdata : str
        The input string to be converted.
    threads : int
        The number of threads to use for parallel processing.

    Returns
    -------
    list
        A list of integers representing the ASCII values of the characters in the input string minus 33.

    Examples
    --------
    >>> _qual_to_ord_dispatcher("Hello", 2)
    [32, 29, 34, 34, 37]
    """

    def _create_chunks(lst: str, n: int) -> Generator:
        """
        Splits a list into approximately equal-sized chunks.

        Parameters
        ----------
        lst : list
            The list to be split.
        n : int
            The number of chunks to create.

        Yields
        ------
        generator
            A generator that yields the chunks.

        Examples
        --------
        >>> lst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        >>> for chunk in chunks(lst, 3):
        ...     print(chunk)
        [1, 2, 3, 4]
        [5, 6, 7]
        [8, 9, 10]
        """
        chunk_size, remainder = divmod(len(lst), n)
        start = 0
        for i in range(n):
            end = start + chunk_size
            if i < remainder:
                end += 1
            yield lst[start:end]
            start = end

    qdata_chunks: list[str] = list(_create_chunks(qdata, threads))
    ordinal_quality_list: list[int] = []
    with ProcessPoolExecutor(max_workers=threads) as pool:
        results = pool.map(_process_chunk, qdata_chunks)
        for result in results:
            ordinal_quality_list.extend(result)
    return ordinal_quality_list


def _process_chunk(chunk: str) -> list[int]:
    """
    Converts each character in the given chunk to its corresponding ASCII value minus 33.

    Parameters
    ----------
    chunk : str
        The input chunk of characters.

    Returns
    -------
    list of int
        A list of integers representing the ASCII values of the characters in the chunk minus 33.

    Examples
    --------
    >>> _process_chunk("Hello")
    [32, 29, 34, 34, 37]

    Notes
    -----
    This function uses the `ord` function to convert each character to its corresponding ASCII value.
    The ASCII value is then subtracted by 33 to obtain the desired result.
    """
    return [ord(character) - 33 for character in chunk]


def _determine_preset(
    avg_len: float, avg_qual: float, quality_range: int, length_range: int
) -> str:
    def _determine_sequence_variance(
        quality_range: int, length_range: int, avg_qual: float, avg_len: float
    ) -> Tuple[float, float]:
        """
        Calculate the variability of a sequence based on the difference between the average quality score and the quality score of the
        current sequence, and the difference between the average length and the length of the current
        sequence.

        Parameters
        ----------
        quality_range : int
            The difference between the highest and lowest quality scores.
        length_range : int
            The difference between the longest and shortest sequence.
        avg_qual : float
            The average quality score of the sequence.
        avg_len : float
            The average length of the sequences.

        Returns
        -------
        Tuple[float, float]
            The quality variance and length variance.

        Examples
        --------
        >>> _determine_sequence_variance(20, 10, 30.0, 5.0)
        (0.6666666666666666, 2.0)

        Notes
        -----
        The quality variance is calculated as the ratio of the quality range to the average quality score.
        The length variance is calculated as the ratio of the length range to the average length.

        """
        quality_variancy = quality_range / avg_qual
        length_variancy = length_range / avg_len
        return quality_variancy, length_variancy

    def _determine_sequence_stability(
        quality_range: int, length_range: int, avg_qual: float, avg_len: float
    ) -> float:
        """
        Calculate the stability of a sequence based on the variability of the quality scores and the lengths of the sequences.

        Parameters
        ----------
        quality_range : int
            The difference between the highest and lowest quality scores.
        length_range : int
            The difference between the longest and shortest sequence.
        avg_qual : float
            The average quality score of the sequence.
        avg_len : float
            The average length of the sequences.

        Returns
        -------
        float
            The percentage of stability of the sequence.

        Examples
        --------
        >>> _determine_sequence_stability(20, 10, 30.0, 5.0)
        32.33333333333333

        Notes
        -----
        The stability of a sequence is calculated by subtracting the sum of the quality score variance and length variance from 100.
        """
        quality_variance, length_variance = _determine_sequence_variance(
            quality_range, length_range, avg_qual, avg_len
        )
        return 100 - (quality_variance + length_variance)

    def _is_long_read(avg_len: float) -> bool:
        """
        Determine if a read is considered long based on its average length.

        Parameters
        ----------
        avg_len : float
            The average length of the read.

        Returns
        -------
        bool
            True if the read is considered long (average length > 300), False otherwise.

        Examples
        --------
        >>> _is_long_read(250)
        False
        >>> _is_long_read(350)
        True

        Notes
        -----
        This function considers a read to be long if its average length is greater than 300.
        """
        return avg_len > 300

    if (
        _determine_sequence_stability(quality_range, length_range, avg_qual, avg_len)
        > 97
    ):
        if _is_long_read(avg_len) is False:
            # this is probably 'short read' illumina NextSeq data
            # --> set the 'SR' preset
            return "sr"
        # ! previous if-statement is not False.
        # this is probably 'long read' illumina MiSeq data
        # --> the 'SR' preset still applies but we keep it split
        # in case a custom set of parameters is necessary in the future
        return "sr"
    if _is_long_read(avg_len) is True:
        # this is probably oxford nanopore data
        # --> set the preset to 'map-ont'
        return "map-ont"
    # ! previous if-statement is not True.
    # this might be very 'unstable' nextseq data,
    # or from a platform we currently dont really support officially.
    # fallback to 'sr' preset
    return "sr"


def _calc_avg_read_length(sequence_list: List[str]) -> float:
    """
    Calculate the average length of a list of sequences.

    Parameters
    ----------
    sequence_list : List[str]
        A list of sequences.

    Returns
    -------
    float
        The average length of the sequences in the list.

    Examples
    --------
    >>> _calc_avg_read_length(['ATCG', 'ATCGATCG', 'AT'])
    5.0

    Notes
    -----
    This function calculates the average length of the sequences in the given list.
    The average length is calculated by summing the lengths of all sequences and dividing it by the total number of sequences.

    The function assumes that all sequences in the list are of type str.

    """
    return sum(map(len, sequence_list)) / float(len(sequence_list))


def _calc_avg_read_qual(quality_list: List[int]) -> float:
    """
    Calculate the average quality score of a list of quality scores.

    Parameters
    ----------
    quality_list : List[int]
        A list of quality scores.

    Returns
    -------
    float
        The average quality score of the quality scores in the list.

    Examples
    --------
    >>> _calc_avg_read_qual([20, 30, 40])
    30.0

    Notes
    -----
    This function calculates the average quality score by summing up all the quality scores in the list
    and dividing it by the total number of quality scores.

    The average quality score is a measure of the overall quality of the sequence data.
    Higher average quality scores indicate better quality data.

    """
    return sum(quality_list) / len(quality_list)


def _get_unique_quality_scores(quality_list: List[int]) -> int:
    """
    Return the number of unique quality scores in a list.

    Parameters
    ----------
    quality_list : List[int]
        A list of quality scores.

    Returns
    -------
    int
        The number of unique quality scores in the list.

    Examples
    --------
    >>> GetQualRange([20, 30, 20, 40])
    3
    """
    return len(set(quality_list))


def _get_unique_read_lengths(sequence_list: List[str]) -> int:
    """
    Return the number of unique lengths of strings in a list.

    Parameters
    ----------
    sequence_list : List[str]
        A list of strings.

    Returns
    -------
    int
        The number of unique lengths in the list of strings.

    Examples
    --------
    >>> _get_unique_read_lengths(['ATCG', 'ATCGATCG', 'AT'])
    2

    Notes
    -----
    This function calculates the number of unique lengths of strings in the given list.
    It does this by iterating over each string in the list and getting its length.
    The lengths are then converted into a set to remove duplicates, and the length of the set is returned.
    """
    all_read_lengths = [len(i) for i in sequence_list]
    return len(set(all_read_lengths))


def _sequence_statistics_dispatcher(
    reads_list: List[str], ordinal_qualities_list: List[int], threads: int
) -> Tuple[float, float, int, int]:
    """
    Calculate sequence statistics using multiple threads.

    Parameters
    ----------
    reads_list : List[str]
        A list of DNA sequences.
    ordinal_qualities_list : List[int]
        A list of quality scores corresponding to the DNA sequences.
    threads : int
        The number of threads to use for parallel execution.

    Returns
    -------
    Tuple[float, float, int, int]
        A tuple containing the average length of the DNA sequences,
        the average quality score, the number of unique quality scores, and the number of unique
        DNA sequence lengths.
    """

    with ProcessPoolExecutor(max_workers=threads) as pool:
        future_averagelength = pool.submit(_calc_avg_read_length, reads_list)
        future_averagequal = pool.submit(_calc_avg_read_qual, ordinal_qualities_list)
        future_qualityrange = pool.submit(
            _get_unique_quality_scores, ordinal_qualities_list
        )
        future_lengthrange = pool.submit(_get_unique_read_lengths, reads_list)

        avg_len = future_averagelength.result()
        avg_qual = future_averagequal.result()
        quality_range = future_qualityrange.result()
        length_range = future_lengthrange.result()
    return avg_len, avg_qual, quality_range, length_range
