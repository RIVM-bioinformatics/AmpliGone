from concurrent.futures import ThreadPoolExecutor


def Calculate_avg_seq_len(Slist):
    """Takes a list of sequences, and returns the average length of the sequences

    Parameters
    ----------
    Slist
        a list of sequences

    Returns
    -------
        The average length of the sequences in the list.

    """
    return sum(map(len, Slist)) / float(len(Slist))


def Calculate_avg_seq_qual(Qlist):
    """This function takes a list of quality scores and returns the average quality score

    Parameters
    ----------
    Qlist
        a list of quality scores

    Returns
    -------
        The average quality score of the sequence.

    """
    return sum(Qlist) / len(Qlist)


def GetQualRange(Qlist):
    """Takes a list of quality scores, converts them to a set, and returns the length of that set

    Parameters
    ----------
    Qlist
        a list of quality scores

    Returns
    -------
        The number of unique values in the list.

    """
    return len(list(set(Qlist)))


def GetLenRange(SList):
    """Takes a list of strings, and returns the number of different lengths of strings in the list

    Parameters
    ----------
    SList
        A list of strings.

    Returns
    -------
        The number of unique lengths in the list of strings.

    """
    AllLengths = [len(i) for i in SList]
    return len(list(set(AllLengths)))


def SequenceVariability(diffQuals, difflen, avg_qual, avg_len):
    """The function takes the difference between the average quality score and the quality score of the
    current sequence, and the difference between the average length and the length of the current
    sequence, and returns the ratio of these differences to the average quality score and average
    length, respectively.

    Parameters
    ----------
    diffQuals
        the difference between the highest and lowest quality scores
    difflen
        the difference between the longest and shortest sequence
    avg_qual
        average quality score of the sequence
    avg_len
        average length of the sequences

    Returns
    -------
        the quality variance and length variance.

    """
    QualVariance = diffQuals / avg_qual
    LengthVariance = difflen / avg_len
    return QualVariance, LengthVariance


def SequenceStability(QualRange, LengthRange, avg_qual, avg_len):
    """Function takes in the range of quality scores and the range of sequence lengths, and returns the
    stability of the sequences as a percentage.

    Parameters
    ----------
    QualRange
        The range of quality scores in the sequence.
    LengthRange
        The range of lengths of the sequences in the file.
    avg_qual
        average quality score of the sequence
    avg_len
        average length of the sequences

    Returns
    -------
        The percentage of stability of the sequence.

    """
    QualVar, LenVar = SequenceVariability(QualRange, LengthRange, avg_qual, avg_len)
    return 100 - (QualVar + LenVar)


def IsLongRead(avg_len):
    """If the average read length is greater than 300, then the read is considered to be long

    Parameters
    ----------
    avg_len
        average read length

    Returns
    -------
        A boolean value

    """
    return avg_len > 300


def FindPreset(threads, data):
    """Takes a dataframe with the sequence and quality data, and returns a preset name and a list of
    scoring parameters

    Parameters
    ----------
    threads
        number of threads to use for the calculations
    data
        a pandas dataframe containing the reads and qualities

    Returns
    -------
        a tuple of two values.
    The first value is a string that represents the preset.
    The second value is a list of integers that represent the scoring matrix.

    """

    ReadList = data["Sequence"].tolist()
    QualList = [
        ord(character) - 33
        for character in [
            x for y in [list(item) for item in data["Qualities"].tolist()] for x in y
        ]
    ]

    if len(ReadList) < 1:
        return None

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
