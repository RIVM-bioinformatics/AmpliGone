from functools import lru_cache
from typing import List


@lru_cache(maxsize=2000000)
def PositionInOrBeforePrimer(pos: int, clist: List[int], max_lookaround: int) -> bool:
    """
    Determine if a position is within the maximum distance of the closest position in the list of primer positions
    and the position is less than or equal to the closest position in the list.

    Parameters
    ----------
    pos : int
        The start or stop position of the read.
    clist : List[int]
        A list of primer positions.
    max_lookaround : int
        The maximum distance from the position to the closest primer.

    Returns
    -------
    bool
        True if the position is within the maximum distance of the closest position in the list AND the position is
        less than or equal to the closest position in the list, False otherwise.

    """
    d = lambda x: abs(x - pos)
    near = min(clist, key=d, default=0)
    return abs(pos - near) < max_lookaround and pos <= near


@lru_cache(maxsize=2000000)
def PositionInOrAfterPrimer(pos: int, clist: List[int], max_lookaround: int) -> bool:
    """
    Determine if a position is within the maximum distance of the closest position in the list of primer positions
    and the position is greater than or equal to the closest position in the list.

    Parameters
    ----------
    pos : int
        The start or stop position of the read.
    clist : List[int]
        A list of primer positions.
    max_lookaround : int
        The maximum distance from the position to the closest primer.

    Returns
    -------
    bool
        True if the position is within the maximum distance of the closest position in the list AND the position is
        greater than or equal to the closest position in the list, False otherwise.

    """
    d = lambda x: abs(x - pos)
    near = min(clist, key=d, default=0)
    return abs(pos - near) < max_lookaround and pos >= near
