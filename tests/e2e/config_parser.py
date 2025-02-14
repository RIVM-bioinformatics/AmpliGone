"""
This module provides a configuration parser for reading and parsing YAML configuration files.

Classes
-------
ConfigParser
    A class to parse YAML configuration files and return the configuration as a dictionary.
"""

import yaml

from AmpliGone.log import log as logger


class ConfigParser:
    """
    A class to parse YAML configuration files.

    Attributes
    ----------
    logger : logging.Logger
        The logger instance for logging messages.

    Methods
    -------
    parse_config(config_file: str) -> dict[str, dict[str, dict[str, str]]]
        Parses the given YAML configuration file and returns the configuration as a dictionary.
    """

    def __init__(self) -> None:
        """
        Initializes the ConfigParser with a logger instance.
        """
        self.logger = logger

    def parse_config(self, config_file: str) -> dict[str, dict[str, dict[str, str]]]:
        """
        Parses the given YAML configuration file and returns the configuration as a dictionary.

        Parameters
        ----------
        config_file : str
            The path to the YAML configuration file.

        Returns
        -------
        dict[str, dict[str, dict[str, str]]]
            The parsed configuration as a nested dictionary.
        """
        with open(config_file, "r", encoding="utf-8") as file:
            config: dict[str, dict[str, dict[str, str]]] = yaml.safe_load(file)
            return config

    def edit_config(self, config_file: str) -> None:
        """
        Placeholder for editing configuration files
        """
        raise NotImplementedError("Editing configuration files is not yet supported.")
