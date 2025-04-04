"""
Module for end-to-end testing of the AmpliGone CLI program.

This module contains the TestE2e class which provides methods to run end-to-end tests
on the AmpliGone CLI program. The tests are configured using a YAML file and include
various scenarios such as real-world data, synthetic data, and edge cases.
"""

import logging
import os
from typing import Generator

import pytest

from AmpliGone.__main__ import main
from tests.e2e.config_parser import ConfigParser


class TestE2e:  # pylint: disable=too-few-public-methods
    """
    Class for end-to-end testing of the AmpliGone CLI program.

    This class provides methods to run end-to-end tests on the AmpliGone CLI program.
    The tests are configured using a YAML file and include various scenarios such as
    real-world data, synthetic data, and edge cases.

    Attributes
    ----------
    config_parser : ConfigParser
        An instance of the ConfigParser class to parse the configuration file.
    config : dict
        A dictionary containing the parsed configuration data.

    Methods
    -------
    _cleanup()
        Fixture to clean up output files after tests are run.
    _order_fastq_by_name(fastq_lines)
        Orders the lines of a FASTQ file by read name.
    _compare_outputs(output_file, expected_output_file)
        Compares the output file with the expected output file.
    _args_to_list(args)
        Converts a dictionary of arguments to a list.
    test_ampligone(test_case, caplog)
        Runs the AmpliGone CLI program with the specified test case and captures logs.
    """

    config_parser = ConfigParser()
    config: dict[str, dict[str, dict[str, str]]] = config_parser.parse_config(
        "tests/config.yaml"
    )

    @pytest.fixture(
        scope="session", autouse=True
    )  # session scope is used to run the fixture only once
    def _cleanup(self) -> Generator[None, None, None]:
        """
        Fixture to clean up output files after tests are run.

        This fixture runs after all tests in the session have completed and removes
        any output files generated during the tests.

        Yields
        ------
        None
        """
        yield  # this is to wait for the test to finish before cleaning up
        for case in self.config.values():
            output_path = case["pipeline_args"]["--output"]
            if os.path.exists(output_path):
                os.remove(output_path)

    def _order_fastq_by_name(self, fastq_lines: list[str]) -> list[str]:
        """
        Orders the lines of a FASTQ file by read name.
        Used because the order of reads in a FASTQ file is not guaranteed.

        Parameters
        ----------
        fastq_lines : list of str
            The lines of the FASTQ file.

        Returns
        -------
        list of str
            The ordered lines of the FASTQ file.
        """
        chunks = [fastq_lines[i : i + 4] for i in range(0, len(fastq_lines), 4)]
        sorted_chunks = sorted(chunks, key=lambda x: x[0])
        sorted_fastq = [line for chunk in sorted_chunks for line in chunk]
        return sorted_fastq

    def _compare_outputs(self, output_file: str, expected_output_file: str) -> None:
        """
        Compares the output file with the expected output file.
        Specifically, it compares the FASTQ files generated by the AmpliGone CLI program.

        Parameters
        ----------
        output_file : str
            The path to the output file generated by the AmpliGone CLI program.
        expected_output_file : str
            The path to the expected output file.

        Raises
        ------
        AssertionError
            If the output file does not match the expected output file.
        """
        with (
            open(output_file, "r", encoding="utf-8") as output,
            open(expected_output_file, "r", encoding="utf-8") as expected_output,
        ):
            output_lines = output.readlines()
            expected_output_lines = expected_output.readlines()
            # the fastq file is unordered, so we need to sort it by the read name
            output_lines = self._order_fastq_by_name(output_lines)
            expected_output_lines = self._order_fastq_by_name(expected_output_lines)
            for line1, line2 in zip(output_lines, expected_output_lines):
                if line1 != line2:
                    raise AssertionError(f"Output: {line1}, Expected Output: {line2}")

    def _args_to_list(self, args: dict[str, str]) -> list[str]:
        """
        Converts a dictionary of arguments to a list.

        Parameters
        ----------
        args : dict of str
            The dictionary of arguments.

        Returns
        -------
        list of str
            The list of arguments.
        """
        return [item for pain in args.items() for item in pain if item]

    @pytest.mark.parametrize(
        "test_case", list(config.values()), ids=list(config.keys())
    )
    def test_ampligone(
        self, test_case: dict[str, dict[str, str]], caplog: pytest.LogCaptureFixture
    ) -> None:
        """
        Runs the AmpliGone CLI program with the specified test cases from the config.yml file.


        Parameters
        ----------
        test_case : dict
            The test case configuration.
        caplog : pytest.LogCaptureFixture
            The log capture fixture.

        Raises
        ------
        SystemExit
            If the test case is expected to fail.
        AssertionError
            If the output file does not match the expected output file or if the
            expected log message is not found in the logs.
        """
        args = self._args_to_list(test_case["pipeline_args"])

        with caplog.at_level(logging.DEBUG, logger="rich"):
            if test_case["test_args"]["fails"]:
                with pytest.raises(SystemExit):
                    main(args)
            else:
                main(args)
                self._compare_outputs(
                    test_case["pipeline_args"]["--output"],
                    test_case["test_args"]["comparison_file"],
                )
        if test_case["test_args"].get(
            "expected_log_message", None
        ):  # I dont think its necessary to check the logs in every case
            assert (
                test_case["test_args"]["expected_log_message"].lower()
                in caplog.text.lower()
            )
