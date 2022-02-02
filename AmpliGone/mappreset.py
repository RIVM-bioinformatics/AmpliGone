from concurrent.futures import ThreadPoolExecutor


def Calculate_avg_seq_len(Slist):
    return sum(map(len, Slist)) / float(len(Slist))


def Calculate_avg_seq_qual(Qlist):
    return sum(Qlist) / len(Qlist)


def GetQualRange(Qlist):
    return len(list(set(Qlist)))


def GetLenRange(SList):
    AllLengths = []
    for i in SList:
        AllLengths.append(len(i))
    return len(list(set(AllLengths)))


def SequenceVariability(diffQuals, difflen, avg_qual, avg_len):
    QualVariance = diffQuals / avg_qual
    LengthVariance = difflen / avg_len
    return QualVariance, LengthVariance


def SequenceStability(QualRange, LengthRange, avg_qual, avg_len):
    QualVar, LenVar = SequenceVariability(QualRange, LengthRange, avg_qual, avg_len)
    return 100 - (QualVar + LenVar)


def IsLongRead(avg_len):
    if avg_len > 300:
        return True
    return False


def FindPreset(threads, data):

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
