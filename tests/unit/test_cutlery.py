"""
Unit tests for the `cutlery` module in the AmpliGone package.

This module contains unit tests for the functions in the `cutlery` module, which provide utility functions for processing sequencing reads.

Test Scenarios
--------------
- Test the `position_in_or_before_primer` function with various read positions.
- Test the `position_in_or_after_primer` function with various read positions.

Functions
---------
test_position_in_or_before_primer(read, result)
    Tests the `position_in_or_before_primer` function with various read positions and expected results.
test_position_in_or_after_primer(read, result)
    Tests the `position_in_or_after_primer` function with various read positions and expected results.
"""

import pytest

from AmpliGone import cutlery


@pytest.mark.parametrize(
    "read, result",
    [
        (22, True),
        (25, True),
        (10, False),
        (26, False),
    ],
    ids=[
        "before_within_lookaround",
        "same_as_primerstart",
        "before_outside_lookaround",
        "higher_than_primerstart",
    ],
)
def test_position_in_or_before_primer(read: int, result: bool) -> None:
    """
    Tests the `position_in_or_before_primer` function with various read positions and expected results.

    Parameters
    ----------
    read : int
        The read position to test.
    result : bool
        The expected result indicating whether the read position is in or before the primer.

    Returns
    -------
    None
        This function does not return any value. It asserts that the function's output matches the expected result.
    """
    primer_positions = (25, 35)
    max_lookaround = 10
    outcome = cutlery.position_in_or_before_primer(
        read, primer_positions, max_lookaround
    )
    if outcome != result:
        raise AssertionError(
            f"Expected {result} but got {outcome} while running cutlery.position_in_or_before_primer({read}, {primer_positions}, {max_lookaround})"
        )


@pytest.mark.parametrize(
    "read, result",
    [
        (22, False),
        (10, False),
        (25, True),
        (26, True),
    ],
    ids=[
        "before_within_lookaround",
        "before_outside_lookaround",
        "same_as_primerstart",
        "higher_than_primerstart",
    ],
)
def test_postition_in_or_after_primer(read: int, result: bool) -> None:
    """
    Tests the `position_in_or_after_primer` function with various read positions and expected results.

    Parameters
    ----------
    read : int
        The read position to test.
    result : bool
        The expected result indicating whether the read position is in or after the primer.

    Returns
    -------
    None
        This function does not return any value. It asserts that the function's output matches the expected result.
    """
    primer_positions = (25, 35)
    max_lookaround = 10
    outcome = cutlery.position_in_or_after_primer(
        read, primer_positions, max_lookaround
    )
    if outcome != result:
        raise AssertionError(
            f"Expected {result} but got {outcome} while running cutlery.position_in_or_before_primer({read}, {primer_positions}, {max_lookaround})"
        )
