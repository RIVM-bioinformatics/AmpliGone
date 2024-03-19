import sys
from concurrent.futures import ThreadPoolExecutor
from typing import List, Set, Tuple

import pandas as pd

from AmpliGone.func import log


def Calculate_avg_seq_len(sequence_list: List[str]) -> float:
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
    >>> Calculate_avg_seq_len(['ATCG', 'ATCGATCG', 'AT'])
    5.0
    """
    return sum(map(len, sequence_list)) / float(len(sequence_list))


def Calculate_avg_seq_qual(quality_list: List[int]) -> float:
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
    >>> Calculate_avg_seq_qual([20, 30, 40])
    30.0
    """
    return sum(quality_list) / len(quality_list)


def GetQualRange(quality_list: List[int]) -> int:
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


def GetLenRange(sequence_list: List[str]) -> int:
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
    >>> GetLenRange(['ATCG', 'ATCGATCG', 'AT'])
    2
    """
    AllLengths = [len(i) for i in sequence_list]
    return len(set(AllLengths))


def SequenceVariability(
    QualRange: int, LengthRange: int, avg_qual: float, avg_len: float
) -> Tuple[float, float]:
    """
    Calculate the variability of a sequence based on the difference between the average quality score and the quality score of the
    current sequence, and the difference between the average length and the length of the current
    sequence.

    Parameters
    ----------
    QualRange : int
        The difference between the highest and lowest quality scores.
    LengthRange : int
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
    >>> SequenceVariability(20, 10, 30.0, 5.0)
    (0.6666666666666666, 2.0)
    """
    QualVariance = QualRange / avg_qual
    LengthVariance = LengthRange / avg_len
    return QualVariance, LengthVariance


def SequenceStability(
    QualRange: int, LengthRange: int, avg_qual: float, avg_len: float
) -> float:
    """
    Calculate the stability of a sequence based on the variability of the quality scores and the lengths of the sequences.

    Parameters
    ----------
    QualRange : int
        The difference between the highest and lowest quality scores.
    LengthRange : int
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
    >>> SequenceStability(20, 10, 30.0, 5.0)
    32.33333333333333
    """
    QualVar, LenVar = SequenceVariability(QualRange, LengthRange, avg_qual, avg_len)
    return 100 - (QualVar + LenVar)


def IsLongRead(avg_len: float) -> bool:
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
    >>> IsLongRead(250)
    False
    >>> IsLongRead(350)
    True
    """
    return avg_len > 300


def FindPreset(threads: int, data: pd.DataFrame) -> str:
    """
    Takes a dataframe with the sequence and quality data, and returns a preset name and a list of scoring parameters.

    Parameters
    ----------
    threads : int
        Number of threads to use for the calculations.
    data : pandas.DataFrame
        A dataframe containing the reads and qualities.

    Returns
    -------
    str
        A string that represents the preset.

    Notes
    -----
    This function takes a dataframe with the sequence and quality data, and returns a preset name and a list of scoring parameters.
    It first extracts the sequence and quality data from the dataframe, and then calculates the average length and quality of the sequences,
    as well as the range of quality scores and sequence lengths. Based on these values, it determines the appropriate preset to use for the data.
    If the data is considered stable and short, it returns the 'sr' preset. If the data is considered stable and long, it returns the 'map-ont' preset.
    If the data is considered unstable or from an unsupported platform, it falls back to the 'sr' preset.

    Examples
    --------
    >>> data = pd.DataFrame({'Sequence': ['ATCG', 'GCTA'], 'Qualities': ['!@#$', 'abcd']})
    >>> FindPreset(threads=2, data=data)
    'sr'
    """
    ReadList: List[str] = data["Sequence"].tolist()
    QualList: List[int] = [
        ord(character) - 33
        for character in [
            x for y in [list(item) for item in data["Qualities"].tolist()] for x in y
        ]
    ]

    with ThreadPoolExecutor(max_workers=threads) as ex:
        TP_averagelength = ex.submit(Calculate_avg_seq_len, ReadList)
        TP_averagequal = ex.submit(Calculate_avg_seq_qual, QualList)
        TP_qualityrange = ex.submit(GetQualRange, QualList)
        TP_lengthrange = ex.submit(GetLenRange, ReadList)

        avg_len = TP_averagelength.result()
        avg_qual = TP_averagequal.result()
        QualRange = TP_qualityrange.result()
        LengthRange = TP_lengthrange.result()
    if SequenceStability(QualRange, LengthRange, avg_qual, avg_len) > 97:
        if IsLongRead(avg_len) is False:
            # this is probably 'short read' illumina NextSeq data
            # --> set the 'SR' preset
            return "sr"
        ##! previous if-statement is not False.
        # this is probably 'long read' illumina MiSeq data
        # --> the 'SR' preset still applies but we keep it split
        # in case a custom set of parameters is necessary in the future
        return "sr"
    if IsLongRead(avg_len) is True:
        # this is probably oxford nanopore data
        # --> set the preset to 'map-ont'
        return "map-ont"
    ##! previous if-statement is not True.
    # this might be very 'unstable' nextseq data,
    # or from a platform we currently dont really support officially.
    # fallback to 'sr' preset
    return "sr"


def valid_scoring_list_length(input_list: List[str]) -> bool:
    """
    Check if the length of the input list is either 4, 6, or 7.

    Parameters
    ----------
    input_list : list of str
        The list to check the length of.

    Returns
    -------
    bool
        True if the length of the list is 4, 6, or 7, False otherwise.
    """
    return len(input_list) in {4, 6, 7}


def valid_scoring_elements(input_list: Set[str], required_elements: Set[str]) -> bool:
    """
    Check if all elements of the input list are present in the set of required elements.

    Parameters
    ----------
    input_list : set of str
        The set to check if its elements are in the required elements set.
    required_elements : set of str
        The set of required elements.

    Returns
    -------
    bool
        True if all elements of the input list are in the required elements set, False otherwise.
    """
    return input_list.issubset(required_elements)


def scoring_has_negative_values(input_list: List[int]) -> bool:
    """
    Check if the input list contains any negative values.

    Parameters
    ----------
    input_list : list of int
        The list to check for negative values.

    Returns
    -------
    bool
        True if any item in the list is less than 0, False otherwise.
    """
    return any(item < 0 for item in input_list)


def parse_scoring_matrix(input_matrix: List[str]) -> List[int]:
    """
    Parse the scoring matrix from a list of 'key=value' strings to a list of integers matching the required scoring matrix order.

    Parameters
    ----------
    input_matrix : list of str
        The scoring matrix as a list of strings, where each string is a key-value pair separated by '='.

    Returns
    -------
    list of int
        The scoring matrix as a list of integers, ordered to fit the mappy input.

    Raises
    ------
    SystemExit
        If the input matrix is invalid (e.g., wrong length, contains negative values, invalid keys).

    """
    required_4 = sorted(["match", "mismatch", "gap_o1", "gap_e1"])
    required_6 = sorted(required_4 + ["gap_o2", "gap_e2"])
    required_7 = sorted(required_6 + ["mma"])

    matrix_dict = {
        item.split("=")[0]: int(item.split("=")[1]) for item in sorted(input_matrix)
    }

    if valid_scoring_list_length(list(matrix_dict.keys())) is False:
        log.error(
            f"Invalid scoring matrix length. The scoring-matrix must have a length of 4, 6 or 7 parameters. \nThe following input parameters were given: '[red]{' '.join(input_matrix)}[/red]'. \nAfter parsing, these inputs result in the following: [red]{matrix_dict}[/red]. \nPlease note that adding the same key multiple times will result in the last value being used."
        )
        sys.exit(1)

    # check if the values of the scoring matrix are non-negative integers
    if scoring_has_negative_values(list(matrix_dict.values())) is True:
        log.error(
            "Given scoring matrix contains a negative value. \nThe scoring matrix may only contain non-negative integers. Please check your input and try again."
        )
        sys.exit(1)

    # this section is quite redundant and the same thing is being done multiple times, will see to optimize it later but this is reasonably valid for now.
    ordered_vals: List[int] = []
    matrix_keys = list(matrix_dict.keys())
    if len(matrix_keys) == 4:
        if not set(matrix_keys).issubset(required_4):
            log.error(
                f"Invalid combination of scoring matrix keys. A total of 4 valid scoring-matrix keys were given. \nThe following keys are supported for 4 scoring-matrix keys: '[green]{' | '.join(required_4)}[/green]'. \nThe following keys were given: '[red]{' | '.join(matrix_keys)}[/red]'."
            )
            sys.exit(1)
        ordered_vals = [matrix_dict[key] for key in required_4]
    elif len(matrix_keys) == 6:
        if not set(matrix_keys).issubset(required_6):
            log.error(
                f"Invalid combination of scoring matrix keys. A total of 6 valid scoring-matrix keys were given. \nThe following keys are supported for 6 scoring-matrix keys: '[green]{' | '.join(required_6)}[/green]'. \nThe following keys were given: '[red]{' | '.join(matrix_keys)}[/red]'."
            )
            sys.exit(1)
        ordered_vals = [matrix_dict[key] for key in required_6]
    elif len(matrix_keys) == 7:
        if not set(matrix_keys).issubset(required_7):
            log.error(
                f"Invalid combination of scoring matrix keys. A total of 7 valid scoring-matrix keys were given. \nThe following keys are supported for 7 scoring-matrix keys: '[green]{' | '.join(required_7)}[/green]'. \nThe following keys were given: '[red]{' | '.join(matrix_keys)}[/red]'."
            )
            sys.exit(1)
        ordered_vals = [matrix_dict[key] for key in required_7]
    return ordered_vals
