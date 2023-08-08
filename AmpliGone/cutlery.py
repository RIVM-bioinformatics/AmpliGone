from functools import lru_cache


@lru_cache(maxsize=2000000)
def PositionInOrBeforePrimer(pos, clist, max_lookaround):
    """Given a position, a list of positions, and a maximum distance, return True if the position is
    within the maximum distance of the closest position in the list AND the position is less than or
    equal to the closest position in the list

    Parameters
    ----------
    pos
        the position of the variant
    clist
        a list of primer positions
    max_lookaround
        the maximum distance from the position to the closest primer

    Returns
    -------
        A boolean value.

    """
    d = lambda x: abs(x - pos)
    near = min(clist, key=d, default=0)
    return abs(pos - near) < max_lookaround and pos <= near


@lru_cache(maxsize=2000000)
def PositionInOrAfterPrimer(pos, clist, max_lookaround):
    """Given a position, a list of positions, and a maximum distance, return True if the position is
    within the maximum distance of the closest position in the list AND the position is greater than or equal
    to the closest position in the list

    Parameters
    ----------
    pos
        the start or stop position of the read
    clist
        a list of primer positions
    max_lookaround
        the maximum distance from the position to the closest primer

    Returns
    -------
        A boolean value.

    """
    d = lambda x: abs(x - pos)
    near = min(clist, key=d, default=0)
    return abs(pos - near) < max_lookaround and pos >= near


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
