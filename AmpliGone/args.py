"""
This module provides functionality for parsing command-line arguments for the AmpliGone tool using the argparse library.
It includes custom argument validation functions, a flexible argument formatter, and a rich argument parser for enhanced
command-line interface (CLI) experience.

Functions
---------
get_args(givenargs: List[str]) -> argparse.Namespace
    Parses the given command-line arguments and returns them as an argparse namespace.

Classes
-------
FlexibleArgFormatter(argparse.HelpFormatter)
    A subclass of argparse.HelpFormatter that improves the formatting of help text.

RichParser(argparse.ArgumentParser)
    A subclass of argparse.ArgumentParser that uses rich.print for displaying messages.

"""

import argparse
import multiprocessing
import os
import pathlib
import re
import shutil
import textwrap
from typing import IO, Iterable, List, Optional

import rich

from AmpliGone import __prog__, __version__


def get_args(givenargs: List[str]) -> argparse.Namespace:
    """It takes the arguments given to the script and parses them into the argparse namespace

    Parameters
    ----------
    givenargs
        the arguments given to the script

    Returns
    -------
        The arguments that are given to the program.

    """

    def check_file_extensions(allowed_extensions: Iterable[str], fname: str) -> str:
        """
        Check if the file extension of the file name passed to it is one of the extensions in the list passed to it.

        Parameters
        ----------
        allowed_extensions : iterable of str
            A tuple of strings that are valid extensions.
        fname : str
            The name of the file to be read.

        Returns
        -------
        str
            The absolute path of the file passed.

        Raises
        ------
        ArgumentParser.error
            If the file extension of the file name passed to it is not one of the extensions in the list passed to it.
        """
        ext = "".join(pathlib.Path(fname).suffixes)
        if not any(ext.endswith(c) for c in allowed_extensions):
            parser.error(f"File {fname} doesn't end with one of {allowed_extensions}")
        return os.path.abspath(fname)

    def check_file_exists(fname: str) -> str:
        """Check if the given file `fname` exists and return the absolute path.

        Parameters
        ----------
        fname : str
            The name of the file to be checked.

        Returns
        -------
        str
            The absolute path of the file if it exists.

        Raises
        ------
        ArgumentParser.error
            If the file does not exist.

        """
        if os.path.isfile(fname):
            return fname
        raise argparse.ArgumentTypeError(f'File "{fname}" does not exist.')

    parser = RichParser(
        prog=f"[bold]{__prog__}[/bold]",
        usage=f"[bold]{__prog__}[/bold] \\[required options] \\[optional arguments]",
        description=f"[bold underline]{__prog__}[/bold underline]: An accurate and efficient tool to remove primers from NGS reads in reference-based experiments",
        formatter_class=FlexibleArgFormatter,
        add_help=False,
    )

    # set the number of threads to 2 if the system, has 2 or more threads.
    standard_threads = 2 if multiprocessing.cpu_count() >= 2 else None

    required_args = parser.add_argument_group("Required Arguments")

    required_args.add_argument(
        "--input",
        "-i",
        type=lambda s: check_file_exists(
            check_file_extensions((".fastq", ".fq", ".bam", ".fastq.gz", ".fq.gz"), s)
        ),
        metavar="File",
        help="Input file with reads in either FastQ or BAM format.",
        required=True,
    )

    required_args.add_argument(
        "--output",
        "-o",
        type=lambda s: check_file_extensions(
            (".fastq", ".fq", ".fastq.gz", ".fq.gz"), s
        ),
        metavar="File",
        help="Output (FastQ) file with cleaned reads.",
        required=True,
    )

    required_args.add_argument(
        "--reference",
        "-ref",
        type=lambda s: check_file_exists(check_file_extensions((".fasta", ".fa"), s)),
        metavar="File",
        help="Input Reference genome in FASTA format",
        required=True,
    )
    required_args.add_argument(
        "--primers",
        "-pr",
        type=lambda s: check_file_exists(
            check_file_extensions((".fasta", ".fa", ".bed"), s)
        ),
        metavar="File",
        help="""Used primer sequences in FASTA format or primer coordinates in BED format.\n Note that using bed-files overrides error-rate and ambiguity functionality""",
        required=True,
    )

    optional_args = parser.add_argument_group("Optional Arguments")

    optional_args.add_argument(
        "--amplicon-type",
        "-at",
        default="end-to-end",
        choices=("end-to-end", "end-to-mid", "fragmented"),
        help="Define the amplicon-type, either being [green]'end-to-end'[/green], [green]'end-to-mid'[/green], or [green]'fragmented'[/green]. See the docs for more info :book:",
        required=False,
        metavar="'end-to-end'/'end-to-mid'/'fragmented'",
    )

    optional_args.add_argument(
        "--virtual-primers",
        "-vp",
        default=False,
        action="store_true",
        required=False,
        help="If set, primers closely positioned to each other in the same orientation will be virtually combined into a single primer, ensuring proper removal of all primer-related data in a specific region. This is useful for amplicons that share multiple (alternative) primers for increased specificity.",
    )

    optional_args.add_argument(
        "--fragment-lookaround-size",
        "-fls",
        required=False,
        type=int,
        metavar="N",
        help="The number of bases to look around a primer-site to consider it part of a fragment. Only used if amplicon-type is 'fragmented'. Default is 10",
    )

    optional_args.add_argument(
        "--error-rate",
        "-er",
        type=float,
        default=0.1,
        metavar="N",
        help="The maximum allowed error rate (as a percentage) for the primer search. Use 0 for exact primer matches. (0.1 = 10%% error rate)\n Note that this is only used if the primer-file is in FASTA format. If the primer-file is in BED format, this option is ignored.",
        required=False,
    )

    optional_args.add_argument(
        "--alignment-preset",
        "-ap",
        type=str,
        default=None,
        choices=("sr", "map-ont", "map-pb", "splice"),
        help="The preset to use for alignment of reads against the reference. This can be either 'sr', 'map-ont', 'map-pb', or 'splice'. The alignment-preset can be combined with a custom alignment-scoring matrix. See the docs for more info :book:",
        required=False,
        metavar="'sr'/'map-ont'/'map-pb'/'splice'",
    )

    optional_args.add_argument(
        "--alignment-scoring",
        "-as",
        type=str,
        default=None,
        metavar="KEY=VALUE",
        nargs="+",
        help="The scoring matrix to use for alignment of reads. This should be list of key-value pairs, where the key is the scoring-parameter and the value is a positive integer indicating the scoring-value for that parameter. Possible parameters are \n * (1) match\n * (2) mismatch\n * (3) gap_o1\n * (4) gap_e1\n * (5) gap_o2 (Optional: requires 1,2,3,4)\n * (6) gap_e2 (Optional, requires 1,2,3,4,5)\n * (7) mma (Optional, requires 1,2,3,4,5,6)\nFor example:\n --alignment-scoring match=4 mismatch=3 gap_o1=2 gap_e1=1\nSee the docs for more info :book:",
        required=False,
    )

    optional_args.add_argument(
        "--export-primers",
        "-ep",
        type=lambda s: check_file_extensions((".bed",), s),
        metavar="File",
        help="Output BED file with found primer coordinates if they are actually cut from the reads",
        required=False,
    )

    optional_args.add_argument(
        "--threads",
        "-t",
        type=int,
        default=standard_threads,
        metavar="N",
        help=f"""Number of threads you wish to use.\n Default is the number of available threads in your system ({standard_threads})""",
    )

    optional_args.add_argument(
        "-to",
        action="store_true",
        help=f"If set, {__prog__} will always create the output files even if there is nothing to output. (for example when an empty input-file is given)\n This is useful in (automated) pipelines where you want to make sure that the output files are always created.",
        required=False,
    )

    optional_args.add_argument(
        "--verbose",
        "-V",
        action="store_true",
        help="Prints more information, like DEBUG statements, to the terminal",
        required=False,
    )

    optional_args.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Prints less information, like only WARNING and ERROR statements, to the terminal",
        required=False,
    )

    optional_args.add_argument(
        "--version",
        "-v",
        action="version",
        version=__version__,
        help=f"You are using [bold]{__prog__}[/bold] version {__version__}",
    )

    optional_args.add_argument(
        "--help",
        "-h",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )

    return parser.parse_args(givenargs)


class FlexibleArgFormatter(argparse.HelpFormatter):
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

    def __init__(self, prog: str) -> None:
        term_width = shutil.get_terminal_size().columns
        max_help_position = min(max(24, term_width // 2), 80)
        super().__init__(prog, max_help_position=max_help_position)

    # action is actually an argparse._StoreAction object, which is a subclass of argparse.Action
    # _StoreAction violates the Liskov Substitution Principle
    # see: https://mypy.readthedocs.io/en/stable/common_issues.html#incompatible-overrides
    def _get_help_string(self, action: argparse.Action) -> str:
        """ """
        help_text = action.help
        if (
            action.default != argparse.SUPPRESS
            and help_text is not None
            and "default" not in help_text.lower()
            and action.default is not None
        ):
            help_text += f"\n  ([underline]default: {str(action.default)}[/underline])"
        if not help_text:
            raise AssertionError("Help text should always be present")
        return help_text

    # see comment above
    def _format_action_invocation(self, action: argparse.Action) -> str:
        """ """
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ", ".join(action.option_strings) + " " + args_string

    def _split_lines(self, text: str, width: int) -> list[str]:
        return self._para_reformat(text, width)

    def _fill_text(self, text: str, width: int, _: str) -> str:
        lines = self._para_reformat(text, width)
        return "\n".join(lines)

    def _indents(self, line: str) -> tuple[int, int]:
        """Return line indent level and "sub_indent" for bullet list text."""
        matched_line = re.match(r"( *)", line)
        if not matched_line:
            raise AssertionError("Line should always match this regex pattern: ( *)")
        indent = len(matched_line.group(1))
        if list_match := re.match(r"( *)(([*\-+>]+|\w+\)|\w+\.) +)", line):
            sub_indent = indent + len(list_match.group(2))
        else:
            sub_indent = indent

        return (indent, sub_indent)

    def _split_paragraphs(self, text: str) -> list[str]:
        """Split text in to paragraphs of like-indented lines."""

        text = textwrap.dedent(text).strip()
        text = re.sub("\n\n[\n]+", "\n\n", text)

        last_sub_indent = None
        paragraphs: list[str] = []
        for line in text.splitlines():
            (indent, sub_indent) = self._indents(line)
            is_text = re.search(r"[^\s]", line) is not None

            if is_text and indent == sub_indent == last_sub_indent:
                paragraphs[-1] += f" {line}"
            else:
                paragraphs.append(line)

            last_sub_indent = sub_indent if is_text else None
        return paragraphs

    def _para_reformat(self, text: str, width: int) -> list[str]:
        """Reformat text, by paragraph."""

        paragraphs: list[str] = []

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


class RichParser(argparse.ArgumentParser):
    """
    A subclass of `argparse.ArgumentParser` that overrides the `_print_message` method to use
    `rich.print` instead of `print`
    """

    def _print_message(self, message: str, file: Optional[IO[str]] = None) -> None:
        return rich.print(message)
