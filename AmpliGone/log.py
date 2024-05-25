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
