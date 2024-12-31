"""
Unit tests for the `get_args` function in the AmpliGone module.

This module contains unit tests for the `get_args` function, which parses command-line arguments for the AmpliGone pipeline.

Test Scenarios
--------------
- Test with the `--help` flag to ensure it raises `SystemExit`.
- Test with no arguments to ensure it raises `SystemExit`.
- Test with invalid arguments to ensure it raises `SystemExit`.
- Test with necessary arguments to ensure they are parsed correctly.

Classes
-------
TestArgs
    A class containing unit tests for the `get_args` function.
"""

import pytest

from AmpliGone.args import get_args
from tests.e2e.config_parser import ConfigParser


class TestArgs:
    """
    Unit tests for the `get_args` function in the AmpliGone module.

    This class contains unit tests for the `get_args` function,
    which parses command-line arguments for the AmpliGone pipeline.

    Attributes
    ----------
    config_parser : ConfigParser
        An instance of the ConfigParser class to parse the configuration file.
    config : dict
        The parsed configuration dictionary.
    pipeline_args : dict
        The dictionary of pipeline arguments from the configuration.
    happy_arg_list : list
        The list of arguments for the happy flow scenario.

    Methods
    -------
    test_help()
        Test that the `get_args` function raises `SystemExit` when the `--help` flag is provided.
    test_no_args()
        Test that the `get_args` function raises `SystemExit` when no arguments are provided.
    test_invalid_args()
        Test that the `get_args` function raises `SystemExit` when invalid arguments are provided.
    test_necessary_args()
        Test that the `get_args` function correctly parses the necessary arguments.
    """

    config_parser = ConfigParser()
    config = config_parser.parse_config("tests/config.yaml")
    pipeline_args = config["happy_sars_cov_2"]["pipeline_args"]
    happy_arg_list = [item for pain in pipeline_args.items() for item in pain]

    def test_help(self) -> None:
        """
        Test that the `get_args` function raises `SystemExit` when the `--help` flag is provided.

        This test ensures that the `get_args` function correctly handles the `--help` flag by raising `SystemExit`.

        Parameters
        ----------
        self : TestArgs
            The instance of the test class.

        Returns
        -------
        None
            This function does not return any value. It asserts that `SystemExit` is raised.
        """
        with pytest.raises(SystemExit):
            get_args(["--help"])

    def test_no_args(self) -> None:
        """
        Test that the `get_args` function raises `SystemExit` when no arguments are provided.

        This test ensures that the `get_args` function correctly handles the case where no arguments are provided by raising `SystemExit`.

        Parameters
        ----------
        self : TestArgs
            The instance of the test class.

        Returns
        -------
        None
            This function does not return any value. It asserts that `SystemExit` is raised.
        """
        with pytest.raises(SystemExit):
            get_args([])

    def test_invalid_args(self) -> None:
        """
        Test that the `get_args` function raises `SystemExit` when invalid arguments are provided.

        This test ensures that the `get_args` function correctly handles invalid arguments by raising `SystemExit`.

        Parameters
        ----------
        self : TestArgs
            The instance of the test class.

        Returns
        -------
        None
            This function does not return any value. It asserts that `SystemExit` is raised.
        """
        with pytest.raises(SystemExit):
            get_args(["--invalid"])

    def test_necessary_args(self) -> None:
        """
        Test that the `get_args` function correctly parses the necessary arguments.

        This test ensures that the `get_args` function correctly parses the necessary arguments and matches them with the expected values.

        Parameters
        ----------
        self : TestArgs
            The instance of the test class.

        Returns
        -------
        None
            This function does not return any value. It asserts that the parsed arguments match the expected values.
        """
        args = vars(
            get_args(self.happy_arg_list)
        )  # get_args returns a Namespace object, so we need to convert it to a dictionary

        for args_key, args_value in args.items():
            for expected_key, expected_value in self.pipeline_args.items():
                if args_key == expected_key[2:]:  # remove the '--' from the key
                    if args_key == "amplicon_type":
                        assert args_value == expected_value
                    else:  # args_value has absolute path, expected_value has relative path
                        assert expected_value in args_value
