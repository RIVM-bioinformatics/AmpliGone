import pytest

import AmpliGone.cutlery as cutlery


@pytest.mark.parametrize(
    "read, result",
    [
        (22, True),
        (25, True),
        (10, False),
        (26, False),
    ],
)
def test_position_in_or_before_primer(read: int, result: bool) -> None:
    primer_positions = (25, 35)
    max_lookaround = 10
    assert (
        cutlery.position_in_or_before_primer(read, primer_positions, max_lookaround)
        == result
    )


@pytest.mark.parametrize(
    "read, result",
    [
        (22, False),
        (10, False),
        (25, True),
        (26, True),
    ],
)
def test_postition_in_or_after_primer(read: int, result: bool) -> None:
    primer_positions = (25, 35)
    max_lookaround = 10
    assert (
        cutlery.position_in_or_after_primer(read, primer_positions, max_lookaround)
        == result
    )
