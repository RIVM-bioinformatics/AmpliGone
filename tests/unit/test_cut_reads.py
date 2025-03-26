"""
Unit tests for the `cut_reads` function in the AmpliGone module.

This module contains unit tests for the `cut_reads` function,
which processes sequencing reads by cutting them based on primer locations and reference mapping.

Extra tests that could be added if deemed necessary:
- Tests with faulty inputs, like an empty DataFrame, invalid sequences or qualities.
- Tests with different parameters, like different scoring thresholds or presets.


Test Scenarios
--------------
- Test with reads that are too short to be processed.
- Test with no primers available for the reference.
- Test with valid data and different amplicon types.
- Test with primers partially on the read.
- Test with primers not on the read.

Classes
-------
TestCutReads
    A class containing unit tests for the `cut_reads` function.
"""

from collections import defaultdict

import pandas as pd
import pytest

from AmpliGone.cut_reads import cut_reads

# Length of the read should theoretically be ~42bp because the default k-mer size for the short reads preset (-sr) is 21.
# However, in practice this does not work because the -sr preset includes a `-s 40` parameter which means that -
# the scoring threshold is set to 40. This score is very hard to achieve with only 42bp reads.
# This means that for testing purposes reads should be around 50-100bp long.
# To make sure everything works as expected, the read length that will be used for testing is 100bp.

AMPLICON_TYPES = (
    "end-to-end",
    "end-to-mid",
    "fragmented",
)  # cannot be in the class because the pytest.mark.parametrize decorator needs to access it before the class is created


class TestCutReads:
    """
    Unit tests for the `cut_reads` function in the AmpliGone module.

    This class contains unit tests for the `cut_reads` function,
    which processes sequencing reads by cutting them based on primer locations and reference mapping.

    Attributes
    ----------
    HAPPY_SEQ : str
        A sample sequence used for testing.
    HAPPY_QUAL : str
        A sample quality string corresponding to the HAPPY_SEQ.
    reference : str
        The path to the reference genome sequence used for testing.
    preset : str
        The preset used for minimap2 alignment.
    scoring : list of int
        The scoring matrix used for minimap2 alignment.
    fragment_lookaround_size : int
        The number of bases to look around a fragment when cutting reads.

    Methods
    -------
    test_cut_reads_too_short()
        Test that the `cut_reads` function skips reads that are too short to be processed.
    test_cut_reads_no_primers()
        Test that the `cut_reads` function handles cases where no primers are available for the reference.
    test_cut_reads_happy(amplicon_type)
        Test a happy path scenario for the `cut_reads` function with different amplicon types.
    test_cut_reads_primer_half_on_read(amplicon_type)
        Test a scenario where the FW and RV primers are only partially on the read.
    test_cut_reads_wrong_primers(amplicon_type)
        Test a scenario where the primers are not on the read.
    """

    HAPPY_SEQ = "GGAAATTCATTCTAGGGAGTGACGTGGACCCCGGATTGATACAGGATCACATGTAGAAAAGGTAGTCGGACAAGTTACCGCTACCCTCGACCTCGTGGGG"
    HAPPY_QUAL = "EE10%-1#-@7F&@?7(13;-$)A7.7/3(I(.9)&//,$G9?HA'DG=/;3C)2:@C!2/#;8.#7'98AC;FG>E>;E'>'$100G&44763?0,@7I"

    reference = "tests/data/references/synthetic.fasta"
    preset = "sr"
    scoring: list[int] = []
    fragment_lookaround_size = 10000

    def test_cut_reads_too_short(self) -> None:
        """
        Test that the `cut_reads` function skips reads that are too short to be processed.

        This test ensures that reads shorter than the minimum required length (42bp) are not processed
        and are skipped by the `cut_reads` function. It might be useful to add a warning message in the future,
        for reads under 100bp.

        This does not raise an error, but returns an empty DataFrame.

        Parameters
        ----------
        self : TestCutReads
            The instance of the test class containing the test data and parameters.

        Returns
        -------
        None
            This function does not return any value. It asserts that the result DataFrame is empty.
        """

        short_seq = self.HAPPY_SEQ[:20] + self.HAPPY_SEQ[-20:]
        short_qual = self.HAPPY_QUAL[:20] + self.HAPPY_QUAL[-20:]

        data: tuple[pd.DataFrame, int] = (
            pd.DataFrame(
                {
                    "Readname": ["read_number_1"],
                    "Sequence": [short_seq],
                    "Qualities": [short_qual],
                }
            ),
            0,
        )
        primer_1 = set(range(0, 5))
        primer_2 = set(range(15, 20))
        primer_sets = (
            defaultdict(set, {"synthetic_reference": primer_1}),
            defaultdict(set, {"synthetic_reference": primer_2}),
        )
        result = cut_reads(
            data,
            primer_sets,
            self.reference,
            self.preset,
            self.scoring,
            self.fragment_lookaround_size,
            "end-to-end",
        )
        assert result.empty

    def test_cut_reads_no_primers(self) -> None:
        """
        Test that the `cut_reads` function handles cases where no primers are available for the reference.

        This test ensures that the `cut_reads` function can process reads even when no primers are available for the reference.
        It verifies that the function does not remove any coordinates in such cases.

        Parameters
        ----------
        self : TestCutReads
            The instance of the test class containing the test data and parameters.

        Returns
        -------
        None
            This function does not return any value. It asserts that the "Removed_coordinates" column in the result DataFrame is empty.
        """
        data: tuple[pd.DataFrame, int] = (
            pd.DataFrame(
                {
                    "Readname": ["read_number_1"],
                    "Sequence": [self.HAPPY_SEQ],
                    "Qualities": [self.HAPPY_QUAL],
                }
            ),
            0,
        )
        empty_primer_1: set = set()
        empty_primer_2: set = set()
        primer_sets = (
            defaultdict(set, {"synthetic_reference": empty_primer_1}),
            defaultdict(set, {"synthetic_reference": empty_primer_2}),
        )
        result = cut_reads(
            data,
            primer_sets,
            self.reference,
            self.preset,
            self.scoring,
            self.fragment_lookaround_size,
            "end-to-end",
        )
        assert not result["Removed_coordinates"].iloc[0]  # empty list

    @pytest.mark.parametrize("amplicon_type", AMPLICON_TYPES)
    def test_cut_reads_happy(self, amplicon_type: str) -> None:
        """
        Test a happy path scenario for the `cut_reads` function with different amplicon types.

        This test ensures that the `cut_reads` function processes reads correctly for different amplicon types,
        including "end-to-end", "end-to-mid", and "fragmented".

        Parameters
        ----------
        self : TestCutReads
            The instance of the test class containing the test data and parameters.
        amplicon_type : str
            The type of amplicon, either "end-to-end", "end-to-mid", or "fragmented".

        Returns
        -------
        None
            This function does not return any value. It asserts that the removed coordinates match the expected coordinates.
        """
        data: tuple[pd.DataFrame, int] = (
            pd.DataFrame(
                {
                    "Readname": ["read_number_1"],
                    "Sequence": [self.HAPPY_SEQ],
                    "Qualities": [self.HAPPY_QUAL],
                }
            ),
            0,
        )
        primer_1 = set(range(0, 10))
        primer_2 = set(range(91, 101))
        primer_sets = (
            defaultdict(set, {"synthetic_reference": primer_1}),
            defaultdict(set, {"synthetic_reference": primer_2}),
        )
        result = cut_reads(
            data,
            primer_sets,
            self.reference,
            self.preset,
            self.scoring,
            self.fragment_lookaround_size,
            amplicon_type,
        )
        if amplicon_type == "end-to-end" or amplicon_type == "fragmented":
            ete_coords: list[int] = result["Removed_coordinates"].iloc[0]
            ete_expected_coords = list(primer_1) + list(primer_2)
            assert ete_coords.sort() == ete_expected_coords.sort()
        else:
            assert amplicon_type == "end-to-mid"
            etm_coords: list[int] = result["Removed_coordinates"].iloc[0]
            etm_expected_coords = list(primer_1)
            assert etm_coords.sort() == etm_expected_coords.sort()

    @pytest.mark.parametrize("amplicon_type", AMPLICON_TYPES)
    def test_cut_reads_primer_half_on_read(self, amplicon_type: str) -> None:
        """
        Test a scenario where the FW and RV primers are only partially on the read.

        This test ensures that the `cut_reads` function can handle cases where the forward (FW) and reverse (RV) primers
        are only partially present on the read sequence.

        Parameters
        ----------
        self : TestCutReads
            The instance of the test class containing the test data and parameters.
        amplicon_type : str
            The type of amplicon, either "end-to-end", "end-to-mid", or "fragmented".

        Returns
        -------
        None
            This function does not return any value. It asserts that the removed coordinates match the expected coordinates.
        """
        seq = "ATCGC" + self.HAPPY_SEQ[5:-5] + "ATCGC"
        data: tuple[pd.DataFrame, int] = (
            pd.DataFrame(
                {
                    "Readname": ["read_number_1"],
                    "Sequence": [seq],
                    "Qualities": [self.HAPPY_QUAL],
                }
            ),
            0,
        )
        primer_1 = set(range(0, 10))
        primer_2 = set(range(91, 101))
        primer_sets = (
            defaultdict(set, {"synthetic_reference": primer_1}),
            defaultdict(set, {"synthetic_reference": primer_2}),
        )
        result = cut_reads(
            data,
            primer_sets,
            self.reference,
            self.preset,
            self.scoring,
            self.fragment_lookaround_size,
            amplicon_type,
        )
        if amplicon_type == "end-to-end" or amplicon_type == "fragmented":
            ete_coords: list[int] = result["Removed_coordinates"].iloc[0]
            ete_expected_coords = list(primer_1) + list(primer_2)
            assert ete_coords.sort() == ete_expected_coords.sort()
        else:
            assert amplicon_type == "end-to-mid"
            etm_coords: list[int] = result["Removed_coordinates"].iloc[0]
            etm_expected_coords = list(primer_1)
            assert etm_coords.sort() == etm_expected_coords.sort()

    @pytest.mark.parametrize("amplicon_type", AMPLICON_TYPES)
    def test_cut_reads_wrong_primers(self, amplicon_type: str) -> None:
        """
        Test a scenario where the primers are not on the read.

        This test ensures that the `cut_reads` function can handle cases where the forward (FW) and reverse (RV) primers
        are not present on the read sequence.

        Parameters
        ----------
        self : TestCutReads
            The instance of the test class containing the test data and parameters.
        amplicon_type : str
            The type of amplicon, either "end-to-end", "end-to-mid", or "fragmented".

        Returns
        -------
        None
            This function does not return any value. It asserts that no coordinates are removed.
        """
        seq = "ATCGCATCGC" + self.HAPPY_SEQ[10:-10] + "ATCGCATCGC"
        data: tuple[pd.DataFrame, int] = (
            pd.DataFrame(
                {
                    "Readname": ["read_number_1"],
                    "Sequence": [seq],
                    "Qualities": [self.HAPPY_QUAL],
                }
            ),
            0,
        )
        primer_1 = set(range(0, 10))
        primer_2 = set(range(91, 101))
        primer_sets = (
            defaultdict(set, {"synthetic_reference": primer_1}),
            defaultdict(set, {"synthetic_reference": primer_2}),
        )
        result = cut_reads(
            data,
            primer_sets,
            self.reference,
            self.preset,
            self.scoring,
            self.fragment_lookaround_size,
            amplicon_type,
        )
        assert not result["Removed_coordinates"].iloc[0]
