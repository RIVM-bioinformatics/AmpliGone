import logging
import shutil
import textwrap
from argparse import SUPPRESS, ArgumentParser, HelpFormatter
from typing import IO, Optional

import regex as _re
import rich
from rich.highlighter import NullHighlighter
from rich.logging import RichHandler

# Central logging object using Rich's logging library
FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET",
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


class FlexibleArgFormatter(HelpFormatter):
    """
    A subclass of ArgumentParser.HelpFormatter that fixes spacing in the help text and respects bullet points.
    Especially useful for multi-line help texts combined with default values.

    This class has taken a lot of inspiration from the 'argparse-formatter' made by Dave Steele; https://github.com/davesteele/argparse_formatter

    Direct link to davesteele's original code that served as an inspiration for this class: https://github.com/davesteele/argparse_formatter/blob/a15d89a99e20b0cad4c389a2aa490c551cef4f9c/argparse_formatter/flexi_formatter.py

    ---
    This class helps to alleviate the following points of the ArgParse help formatting:
    * The help text will be aligned with the argument name/flags, instead of printing the help description on a newline
    * Adjusting the width of the help text in relationship to the width of the terminal to make sure there is enough space between the argument name/flags and the help text (thus not overloading the end-user with an unreadable wall of text)
    * Adding a default value to the help text (on a newline, and indented) if one is provided in the ArgParse constructor
    * Respecting bullet points in the help description
    * Respecting newlines in the help description (you may have to add a space after the newline to make sure it is properly catched by the formatter)
    * Respecting indentation in the help description (up to a certain degree)
    * Changes the behaviour of the metavar to be only printed once per long AND shorthand argument, instead of printing the metavar multiple times for every possible flag.
    """

    def __init__(self, prog):
        term_width = shutil.get_terminal_size().columns
        max_help_position = min(max(24, term_width // 2), 80)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """ """
        help_text = action.help
        if (
            action.default != SUPPRESS
            and "default" not in help_text.lower()
            and action.default is not None
        ):
            help_text += f"\n  ([underline]default: {str(action.default)}[/underline])"
        return help_text

    def _format_action_invocation(self, action):
        """ """
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ", ".join(action.option_strings) + " " + args_string

    def _split_lines(self, text, width):
        return self._para_reformat(text, width)

    def _fill_text(self, text, width, indent):
        lines = self._para_reformat(text, width)
        return "\n".join(lines)

    def _indents(self, line):
        """Return line indent level and "sub_indent" for bullet list text."""

        indent = len(_re.match(r"( *)", line).group(1))
        if list_match := _re.match(r"( *)(([*\-+>]+|\w+\)|\w+\.) +)", line):
            sub_indent = indent + len(list_match.group(2))
        else:
            sub_indent = indent

        return (indent, sub_indent)

    def _split_paragraphs(self, text):
        """Split text in to paragraphs of like-indented lines."""

        text = textwrap.dedent(text).strip()
        text = _re.sub("\n\n[\n]+", "\n\n", text)

        last_sub_indent = None
        paragraphs = []
        for line in text.splitlines():
            (indent, sub_indent) = self._indents(line)
            is_text = _re.search(r"[^\s]", line) is not None

            if is_text and indent == sub_indent == last_sub_indent:
                paragraphs[-1] += f" {line}"
            else:
                paragraphs.append(line)

            last_sub_indent = sub_indent if is_text else None
        return paragraphs

    def _para_reformat(self, text, width):
        """Reformat text, by paragraph."""

        paragraphs = []
        for paragraph in self._split_paragraphs(text):

            (indent, sub_indent) = self._indents(paragraph)

            paragraph = self._whitespace_matcher.sub(" ", paragraph).strip()
            new_paragraphs = textwrap.wrap(
                text=paragraph,
                width=width,
                initial_indent=" " * indent,
                subsequent_indent=" " * sub_indent,
            )

            # Blank lines get eaten by textwrap, put it back with [' ']
            paragraphs.extend(new_paragraphs or [" "])

        return paragraphs


class RichParser(ArgumentParser):
    """
    A subclass of `argparse.ArgumentParser` that overrides the `_print_message` method to use
    `rich.print` instead of `print`
    """

    def _print_message(self, message: str, file: Optional[IO[str]] = None) -> None:
        return rich.print(message)
