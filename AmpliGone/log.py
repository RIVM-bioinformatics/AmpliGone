"""
This module sets up a central logging object using the Rich library for enhanced logging output.

Imports
--------
import logging
    Standard Python logging module.
from rich.highlighter import NullHighlighter
    Import NullHighlighter from the Rich library to disable highlighting.
from rich.logging import RichHandler
    Import RichHandler from the Rich library to handle logging with Rich's features.

Variables
---------
FORMAT : str
    The format string for log messages.
log : logging.Logger
    The central logging object configured to use RichHandler.

Notes
-----
This module configures the logging system to use RichHandler from the Rich library, which provides enhanced logging output with features like rich text formatting and better readability. The logging level is set to DEBUG, and the log messages are formatted to include the date and time.

Examples
--------
>>> from log import log
>>> log.debug("This is a debug message.")
>>> log.info("This is an info message.")
>>> log.warning("This is a warning message.")
>>> log.error("This is an error message.")
>>> log.critical("This is a critical message.")
"""

import logging

from rich.highlighter import NullHighlighter
from rich.logging import RichHandler

# Central logging object using Rich's logging library
FORMAT = "%(message)s"
logging.basicConfig(
    level="INFO",
    format=FORMAT,
    datefmt="[%d/%m/%y %H:%M:%S]",
    handlers=[
        RichHandler(
            show_path=False,
            omit_repeated_times=False,
            markup=True,
            highlighter=NullHighlighter(),
        )
    ],
)
log = logging.getLogger("rich")
