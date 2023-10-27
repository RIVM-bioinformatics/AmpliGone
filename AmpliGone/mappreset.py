from concurrent.futures import ThreadPoolExecutor
from typing import List, Tuple

import pandas as pd


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


def FindPreset(threads: int, data: pd.DataFrame) -> Tuple[str, List[int]]:
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
    Tuple[str, List[int]]
        A tuple of two values. The first value is a string that represents the preset. The second value is a list of integers that represent the scoring matrix.

    Raises
    ------
    None

    Examples
    --------
    >>> data = pd.DataFrame({'Sequence': ['ATCG', 'GCTA'], 'Qualities': ['!@#$', 'abcd']})
    >>> FindPreset(threads=2, data=data)
    ('sr', [])

    Notes
    -----
    This function takes a dataframe with the sequence and quality data, and returns a preset name and a list of scoring parameters.
    It first extracts the sequence and quality data from the dataframe, and then calculates the average length and quality of the sequences,
    as well as the range of quality scores and sequence lengths. Based on these values, it determines the appropriate preset to use for the data.
    If the data is considered stable and short, it returns the 'sr' preset. If the data is considered stable and long, it returns the 'map-ont' preset.
    If the data is considered unstable or from an unsupported platform, it falls back to the 'sr' preset.
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
            return "sr", []
        ##! previous if-statement is not False.
        # this is probably 'long read' illumina MiSeq data
        # --> the 'SR' preset still applies but we keep it split
        # in case a custom set of parameters is necessary in the future
        return "sr", []
    if IsLongRead(avg_len) is True:
        # this is probably oxford nanopore data
        # --> set the preset to 'map-ont'
        O1, O2 = 8, 24
        E1, E2 = 2, 0
        A, B = 4, 4
        scoring = [A, B, O1, E1, O2, E2]
        return "map-ont", scoring
    ##! previous if-statement is not True.
    # this might be very 'unstable' nextseq data,
    # or from a platform we currently dont really support officially.
    # fallback to 'sr' preset
    return "sr", []
