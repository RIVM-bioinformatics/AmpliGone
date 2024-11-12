import os

import pytest

from AmpliGone.args import get_args


@pytest.mark.parametrize(
    "arg_list, expected",
    [
        (
            [
                "--input",
                "tests/data/sars-cov-2/example1/SRR30635841_1.fastq",
                "--primers",
                "tests/data/sars-cov-2/primers/ARTIC-V5.3.2.fasta",
                "--reference",
                "tests/data/sars-cov-2/references/SARS-CoV-2-reference.fasta",
                "--output",
                "tests/data/example1/output.fastq",
            ],
            {
                "input": "tests/data/sars-cov-2/example1/SRR30635841_1.fastq",
                "output": "tests/data/example1/output.fastq",
                "reference": "tests/data/sars-cov-2/references/SARS-CoV-2-reference.fasta",
                "primers": "tests/data/sars-cov-2/primers/ARTIC-V5.3.2.fasta",
                "amplicon_type": "end-to-end",
                "virtual_primers": False,
                "fragment_lookaround_size": None,
                "export_primers": None,
                "threads": 20,
                "to": False,
                "error_rate": 0.1,
                "alignment_preset": None,
                "alignment_scoring": None,
                "verbose": False,
                "quiet": False,
            },
        ),
        (
            ["--help"],
            {
                "result": None,
            },
        ),
    ],
)
def test_get_args(
    arg_list: list[str], expected: dict[str, str | bool | float | int | None]
) -> None:
    if expected.get("result", "not-help") is None:
        # --help will raise SystemExit
        with pytest.raises(SystemExit):
            get_args(arg_list)
        return
    args = get_args(arg_list)
    args_dict = vars(args)

    # vars(args) has absolute paths, so we need to convert the expected paths to absolute
    for key in ["input", "output", "reference", "primers"]:
        path = expected[key]
        if not isinstance(path, str):
            raise AssertionError(f"{key} should be a path string")
        expected[key] = os.path.abspath(path)
    if args_dict != expected:
        raise AssertionError(f"Expected: {expected}, Got: {args_dict}")
