import argparse
import glob
import os
import shutil


class MyHelpFormatter(argparse.RawTextHelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """

    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ["COLUMNS"] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 2), 80)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        help_text = action.help
        if (
            action.default != argparse.SUPPRESS
            and "default" not in help_text.lower()
            and action.default is not None
        ):
            help_text += " (default: " + str(action.default) + ")"
        return help_text


class color:
    PURPLE = "\033[95m"
    CYAN = "\033[96m"
    DARKCYAN = "\033[36m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"
