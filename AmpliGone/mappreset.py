from os import EX_SOFTWARE, readlink, truncate
from .io_ops import ReadsToList
from concurrent.futures import ThreadPoolExecutor


def Calculate_avg_seq_len(Slist):
    return sum(map(len, Slist))/float(len(Slist))

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
    QualVariance = (diffQuals / avg_qual)
    LengthVariance = (difflen / avg_len)
    return QualVariance, LengthVariance

def SequenceStability(QualRange, LengthRange, avg_qual, avg_len):
    QualVar, LenVar = SequenceVariability(QualRange, LengthRange, avg_qual, avg_len)
    return (100 - (QualVar + LenVar))

def IsLongRead(avg_len):
    if avg_len > 300:
        return True
    elif avg_len < 300:
        return False

def FindPreset(inputfile, threads):
    ReadList, QualList = ReadsToList(inputfile)
    
    with ThreadPoolExecutor(max_workers=threads) as exec:
        TP_averagelength = exec.submit(Calculate_avg_seq_len, ReadList)
        TP_averagequal = exec.submit(Calculate_avg_seq_qual, QualList)
        TP_qualityrange = exec.submit(GetQualRange, QualList)
        TP_lengthrange = exec.submit(GetLenRange, ReadList)
        
        avg_len = TP_averagelength.result()
        avg_qual = TP_averagequal.result()
        QualRange = TP_qualityrange.result()
        LengthRange = TP_lengthrange.result()
    
    if SequenceStability(QualRange, LengthRange, avg_qual, avg_len) > 97:
        if IsLongRead(avg_len) is False:
            # this is probably 'short read' illumina NextSeq data --> set the 'SR' preset
            preset = 'sr'
        elif IsLongRead(avg_len) is True:
            # this is probably 'long read' illumina MiSeq data --> the 'SR' preset still applies but we keep it split up in case a custom set of parameters is necessary in the future
            preset = 'sr'
    else:
        if IsLongRead(avg_len) is True:
            # this is probably oxford nanopore data --> set the preset to 'map-ont'
            preset = 'map-ont'
        if IsLongRead(avg_len) is False:
            # this might be very 'unstable' nextseq data, or from a platform we currently dont really support officially. fallback to 'sr' preset
            preset = 'sr'
    return preset
