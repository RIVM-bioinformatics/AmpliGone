from functools import lru_cache


@lru_cache(maxsize=2000000)
def PositionInOrBeforePrimer(pos, clist):
    d = lambda x: abs(x - pos)
    near = min(clist, key=d)
    return pos <= near


@lru_cache(maxsize=2000000)
def PositionInOrAfterPrimer(pos, clist):
    d = lambda x: abs(x - pos)
    near = min(clist, key=d)
    return pos >= near


@lru_cache(maxsize=2000000)
def ReadBeforePrimer(pos, clist):
    if pos in clist:
        return False
    d = lambda x: abs(x - pos)
    near = min(clist, key=d)
    if pos <= near:
        return True
    return False


@lru_cache(maxsize=2000000)
def ReadAfterPrimer(pos, clist):
    if pos in clist:
        return False
    d = lambda x: abs(x - pos)
    near = min(clist, key=d)
    if pos >= near:
        return True
    return False


@lru_cache(maxsize=600000)
def slice_fw_left(start, seq, qual):
    start = start + 1
    tseq = seq[1:]
    tqual = qual[1:]

    return tseq, tqual, start


@lru_cache(maxsize=600000)
def slice_fw_right(end, seq, qual):
    end = end - 1
    tseq = seq[:-1]
    tqual = qual[:-1]

    return tseq, tqual, end


@lru_cache(maxsize=600000)
def slice_rv_left(start, seq, qual):
    start = start + 1
    tseq = seq[:-1]
    tqual = qual[:-1]

    return tseq, tqual, start


@lru_cache(maxsize=600000)
def slice_rv_right(end, seq, qual):
    end = end - 1
    tseq = seq[1:]
    tqual = qual[1:]

    return tseq, tqual, end
