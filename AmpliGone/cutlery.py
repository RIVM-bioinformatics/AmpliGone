"""
This module provides functions to determine if a position is within a specified distance of primer positions.

Functions
---------
position_in_or_before_primer(pos: int, clist: List[int], max_lookaround: int) -> bool
    Determine if a position is within the maximum distance of the closest position in the list of primer positions
    and the position is less than or equal to the closest position in the list.

position_in_or_after_primer(pos: int, clist: List[int], max_lookaround: int) -> bool
    Determine if a position is within the maximum distance of the closest position in the list of primer positions
    and the position is greater than or equal to the closest position in the list.

Notes
-----
These functions use caching to improve performance for repeated calls with the same arguments. The cache size is set to a maximum of 2,000,000 entries.

Examples
--------
>>> from cutlery import position_in_or_before_primer, position_in_or_after_primer
>>> primer_positions = [100, 200, 300]
>>> position_in_or_before_primer(150, primer_positions, 50)
True
>>> position_in_or_after_primer(250, primer_positions, 50)
True
"""

from functools import cache
from typing import List


@cache
def position_in_or_before_primer(
    pos: int, clist: List[int], max_lookaround: int
) -> bool:
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

    def _default(x: int) -> int:
        return abs(x - pos)

    near = min(clist, key=_default, default=0)
    return abs(pos - near) < max_lookaround and pos <= near


@cache
def position_in_or_after_primer(
    pos: int, clist: List[int], max_lookaround: int
) -> bool:
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

    def _default(x: int) -> int:
        return abs(x - pos)

    near = min(clist, key=_default, default=0)
    return abs(pos - near) < max_lookaround and pos >= near
