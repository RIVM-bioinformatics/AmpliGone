"""
This module provides functions to calculate and validate scoring matrices for sequence alignment.

Functions
---------
get_scoring_matrix(input_matrix: Optional[List[str]]) -> List[int]
    Calculate the scoring matrix based on the input matrix.

_input_to_dict(input_matrix: Optional[List[str]]) -> Optional[Dict[str, int]]
    Convert the input matrix to a dictionary.

_valid_scoring_list_length(input_list: List[str]) -> bool
    Check if the length of the input list is either 4, 6, or 7.

_scoring_has_negative_values(input_list: List[int]) -> bool
    Check if the input list contains any negative values.

_sort_matrix_dict(matrix_dict: Dict[str, int], required_4: List[str], required_6: List[str], required_7: List[str]) -> Dict[str, int]
    Sort the given matrix dictionary based on the provided order of keys.

_get_ordered_values(matrix_dict: Dict[str, int], required_4: List[str], required_6: List[str], required_7: List[str]) -> List[int]
    Get the ordered values from the matrix dictionary based on the number of keys present.

_log_invalid_combination_error(matrix_keys: List[str], required_keys: List[str]) -> None
    Log an error message and exit the program when an invalid combination of scoring matrix keys is encountered.

_validate_matrix_combinations(matrix_dict: Dict[str, int]) -> List[int]
    Validate the combinations of matrix values in the given matrix dictionary.

Notes
-----
This module is designed to handle the calculation and validation of scoring matrices used in sequence alignment. It ensures that the input matrices are in the correct format, contain valid values, and have the appropriate length. The module also provides functions to sort and order the matrix values based on predefined requirements.

Examples
--------
>>> input_matrix = ['match=1', 'mismatch=2', 'gap_o1=3', 'gap_e1=4']
>>> get_scoring_matrix(input_matrix)
[1, 2, 3, 4]

>>> input_matrix = ['match=1', 'mismatch=2', 'gap_o1=3', 'gap_e1=4', 'gap_o2=5', 'gap_e2=6']
>>> get_scoring_matrix(input_matrix)
[1, 2, 3, 4, 5, 6]

>>> input_matrix = ['match=1', 'mismatch=2', 'gap_o1=3', 'gap_e1=4', 'gap_o2=5', 'gap_e2=6', 'mma=7']
>>> get_scoring_matrix(input_matrix)
[1, 2, 3, 4, 5, 6, 7]

>>> input_matrix = ['match=1', 'mismatch=2', 'gap_o1=3', 'gap_e1=4', 'gap_o2=5', 'gap_e2=6', 'mma=-7']
>>> get_scoring_matrix(input_matrix)
SystemExit: Given scoring matrix contains a negative value. The scoring matrix may only contain non-negative integers. Please check your input and try again.
"""

import os
import sys
from typing import Dict, List, Optional

from AmpliGone.log import log


def get_scoring_matrix(input_matrix: Optional[List[str]]) -> List[int]:
    """
    Calculate the scoring matrix based on the input matrix.

    Parameters
    ----------
    input_matrix : List[str] or None
        The input matrix used to calculate the scoring matrix. Each element in the list should be in the format 'key=value'.

    Returns
    -------
    List[int]
        The calculated scoring matrix.

    Raises
    ------
    ValueError
        If the input matrix is not in the correct format or contains negative values.

    Notes
    -----
    This function calculates the scoring matrix based on the input matrix.
    The input matrix must have a length of 4, 6, or 7 parameters.
    The scoring matrix may only contain non-negative integers.

    Examples
    --------
    >>> input_matrix = ['match=1', 'mismatch=2', 'gap_o1=3', 'gap_e1=4']
    >>> get_scoring_matrix(input_matrix)
    [1, 2, 3, 4]

    >>> input_matrix = ['match=1', 'mismatch=2', 'gap_o1=3', 'gap_e1=4', 'gap_o2=5', 'gap_e2=6']
    >>> get_scoring_matrix(input_matrix)
    [1, 2, 3, 4, 5, 6]

    >>> input_matrix = ['match=1', 'mismatch=2', 'gap_o1=3', 'gap_e1=4', 'gap_o2=5', 'gap_e2=6', 'mma=7']
    >>> get_scoring_matrix(input_matrix)
    [1, 2, 3, 4, 5, 6, 7]

    >>> input_matrix = ['match=1', 'mismatch=2', 'gap_o1=3', 'gap_e1=4', 'gap_o2=5', 'gap_e2=6', 'mma=7', 'extra=8']
    >>> get_scoring_matrix(input_matrix)
    SystemExit: Invalid scoring matrix length. The scoring-matrix must have a length of 4, 6 or 7 parameters.
    The following input parameters were given: 'match=1 mismatch=2 gap_o1=3 gap_e1=4 gap_o2=5 gap_e2=6 mma=7 extra=8'.
    After parsing, these inputs result in the following: {'match': 1, 'mismatch': 2, 'gap_o1': 3, 'gap_e1': 4, 'gap_o2': 5, 'gap_e2': 6, 'mma': 7, 'extra': 8}.
    Please note that adding the same key multiple times will result in the last value being used.

    >>> input_matrix = ['match=1', 'mismatch=2', 'gap_o1=3', 'gap_e1=4', 'gap_o2=5', 'gap_e2=6', 'mma=-7']
    >>> get_scoring_matrix(input_matrix)
    SystemExit: Given scoring matrix contains a negative value.
    The scoring matrix may only contain non-negative integers. Please check your input and try again.

    """
    log.debug(f"Starting MATRIXPARSER process\t@ ProcessID {os.getpid()}")
    matrix_dict = _input_to_dict(input_matrix)
    if matrix_dict is None:
        return []
    if _valid_scoring_list_length(list(matrix_dict.keys())) is False:
        log.error(
            f"Invalid scoring matrix length. The scoring-matrix must have a length of 4, 6 or 7 parameters. \nThe following input parameters were given: '[red]{' '.join(input_matrix) if input_matrix is not None else None}[/red]'. \nAfter parsing, these inputs result in the following: [red]{matrix_dict}[/red]. \nPlease note that adding the same key multiple times will result in the last value being used."
        )
        sys.exit(1)
    if _scoring_has_negative_values(list(matrix_dict.values())) is True:
        log.error(
            "Given scoring matrix contains a negative value. \nThe scoring matrix may only contain non-negative integers. Please check your input and try again."
        )
        sys.exit(1)
    return _validate_matrix_combinations(matrix_dict)


def _input_to_dict(input_matrix: Optional[List[str]]) -> Optional[Dict[str, int]]:
    """
    Convert the input matrix to a dictionary.

    Parameters
    ----------
    input_matrix : List[str] or None
        The input matrix to be converted. Each element in the list should be in the format 'key=value'.

    Returns
    -------
    dict or None
        A dictionary representation of the input matrix. Returns None if the input matrix is None.

    Raises
    ------
    ValueError
        If the input matrix is not in the correct format.

    Examples
    --------
    >>> input_matrix = ['a=1', 'b=2', 'c=3']
    >>> _input_to_dict(input_matrix)
    {'a': 1, 'b': 2, 'c': 3}

    >>> input_matrix = None
    >>> _input_to_dict(input_matrix)
    None

    >>> input_matrix = ['a=1', 'b=2', 'c']
    >>> _input_to_dict(input_matrix)
    ValueError: Invalid scoring matrix input. The scoring matrix input must be in the format 'key=value'.
            Please check if your input does not contain any spaces and try again.
    """
    if input_matrix is None:
        log.debug("MATRIXPARSER :: No input matrix given.")
        return None
    for item in input_matrix:
        if "=" not in item:
            log.error(
                "Invalid scoring matrix input. The scoring matrix input must be in the format 'key=value'.\nPlease check if your input does not contain any spaces and try again."
            )
            sys.exit(1)
    sort_order = ["match", "mismatch", "gap_o1", "gap_e1", "gap_o2", "gap_e2", "mma"]
    return {
        item.split("=")[0]: int(item.split("=")[1])
        for item in sorted(
            input_matrix, key=lambda x: sort_order.index(x.split("=")[0])
        )
    }


def _valid_scoring_list_length(input_list: List[str]) -> bool:
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

    Raises
    ------
    None

    Examples
    --------
    >>> _valid_scoring_list_length(['A', 'B', 'C', 'D'])
    True

    >>> _valid_scoring_list_length(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
    False
    """
    return len(input_list) in {4, 6, 7}


def _scoring_has_negative_values(input_list: List[int]) -> bool:
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

    Examples
    --------
    >>> _scoring_has_negative_values([1, 2, 3])
    False

    >>> _scoring_has_negative_values([-1, 2, 3])
    True
    """
    return any(item < 0 for item in input_list)


def _sort_matrix_dict(
    matrix_dict: Dict[str, int],
    required_4: List[str],
    required_6: List[str],
    required_7: List[str],
) -> Dict[str, int]:
    """
    Sorts the given matrix dictionary based on the provided order of keys.

    Parameters
    ----------
    matrix_dict : Dict[str, int]
        A dictionary containing key-value pairs, where the keys are strings and the values are integers.
    required_4 : List[str]
        A list of strings representing the required keys that should be placed at the beginning of the sorted dictionary.
    required_6 : List[str]
        A list of strings representing the required keys that should be placed after the required_4 keys in the sorted dictionary.
    required_7 : List[str]
        A list of strings representing the required keys that should be placed after the required_6 keys in the sorted dictionary.

    Returns
    -------
    Dict[str, int]
        A dictionary containing the sorted key-value pairs from the input matrix_dict, based on the provided order of keys.
    """
    order = required_4 + required_6 + required_7
    matrix_dict = {key: matrix_dict[key] for key in order if key in matrix_dict}
    return matrix_dict


def _get_ordered_values(
    matrix_dict: Dict[str, int],
    required_4: List[str],
    required_6: List[str],
    required_7: List[str],
) -> List[int]:
    """
    Get the ordered values from the matrix dictionary based on the number of keys present.

    Parameters
    ----------
    matrix_dict : Dict[str, int]
        A dictionary containing string keys and integer values.
    required_4 : List[str]
        A list of strings representing the required keys when there are 4 keys in the matrix_dict.
    required_6 : List[str]
        A list of strings representing the required keys when there are 6 keys in the matrix_dict.
    required_7 : List[str]
        A list of strings representing the required keys when there are 7 keys in the matrix_dict.

    Returns
    -------
    List[int]
        A list of integers representing the ordered values corresponding to the required keys.

    Raises
    ------
    None

    """
    ordered_vals: List[int] = []
    matrix_keys = list(matrix_dict.keys())
    if len(matrix_keys) == 4:
        if not set(matrix_keys).issubset(required_4):
            _log_invalid_combination_error(matrix_keys, required_4)
        ordered_vals = [matrix_dict[key] for key in required_4]
        log.debug(f"MATRIXPARSER :: Ordered values: {ordered_vals}")
    elif len(matrix_keys) == 6:
        if not set(matrix_keys).issubset(required_6):
            _log_invalid_combination_error(matrix_keys, required_6)
        ordered_vals = [matrix_dict[key] for key in required_6]
        log.debug(f"MATRIXPARSER :: Ordered values: {ordered_vals}")
    elif len(matrix_keys) == 7:
        if not set(matrix_keys).issubset(required_7):
            _log_invalid_combination_error(matrix_keys, required_7)
        ordered_vals = [matrix_dict[key] for key in required_7]
        log.debug(f"MATRIXPARSER :: Ordered values: {ordered_vals}")
    return ordered_vals


def _log_invalid_combination_error(
    matrix_keys: List[str], required_keys: List[str]
) -> None:
    """
    Log an error message and exit the program when an invalid combination of scoring matrix keys is encountered.

    Parameters
    ----------
    matrix_keys : List[str]
        The list of scoring matrix keys that were given.
    required_keys : List[str]
        The list of valid scoring matrix keys that are supported.

    Returns
    -------
    None

    Notes
    -----
    This function logs an error message and exits the program if an invalid combination of scoring matrix keys is encountered.
    The error message includes information about the number of valid scoring matrix keys and the keys that were given.

    Examples
    --------
    >>> matrix_keys = ['A', 'B', 'C']
    >>> required_keys = ['A', 'B', 'C', 'D']
    >>> _log_invalid_combination_error(matrix_keys, required_keys)
    Invalid combination of scoring matrix keys. A total of 4 valid scoring-matrix keys were given.
    The following keys are supported for 4 scoring-matrix keys: 'A | B | C | D'.
    The following keys were given: 'A | B | C'.
    """
    log.error(
        f"Invalid combination of scoring matrix keys. A total of {len(required_keys)} valid scoring-matrix keys were given. \nThe following keys are supported for {len(required_keys)} scoring-matrix keys: '[green]{' | '.join(required_keys)}[/green]'. \nThe following keys were given: '[red]{' | '.join(matrix_keys)}[/red]'."
    )
    sys.exit(1)


def _validate_matrix_combinations(matrix_dict: Dict[str, int]) -> List[int]:
    """
    Validate the combinations of matrix values in the given matrix dictionary.

    Parameters
    ----------
    matrix_dict : Dict[str, int]
        A dictionary containing matrix values.

    Returns
    -------
    List[int]
        A list of ordered matrix values.

    Raises
    ------
    None

    Notes
    -----
    The function validates the combinations of matrix values in the given matrix dictionary.
    It ensures that the required matrix values are present and in the correct order.
    The required matrix values are 'match', 'mismatch', 'gap_o1', 'gap_e1', 'gap_o2', 'gap_e2', and 'mma'.
    The function sorts the matrix dictionary based on the required values and returns the ordered matrix values.

    Examples
    --------
    >>> matrix_dict = {
    ...     'match': 4,
    ...     'mismatch': 4,
    ...     'gap_o1': 6,
    ...     'gap_e1': 2,
    ...     'gap_o2': 14,
    ...     'gap_e2': 1,
    ...     'mma': 10
    ... }
    >>> ordered_values = _validate_matrix_combinations(matrix_dict)
    >>> print(ordered_values)
    [4, 4, 6, 2, 14, 1, 10]
    """
    required_4 = ["match", "mismatch", "gap_o1", "gap_e1"]
    required_6 = required_4 + ["gap_o2", "gap_e2"]
    required_7 = required_6 + ["mma"]

    log.debug(
        "MATRIXPARSER :: Sorting matrix dictionary to fit required alignment matrix."
    )
    matrix_dict = _sort_matrix_dict(matrix_dict, required_4, required_6, required_7)
    log.debug(f"MATRIXPARSER :: {matrix_dict}")
    log.debug("MATRIXPARSER :: Fetching ordered values from sorted matrix dictionary")
    return _get_ordered_values(matrix_dict, required_4, required_6, required_7)
