import os
from typing import Generator

import pytest

from AmpliGone.__main__ import main

INPUT_FILE = "tests/data/sars-cov-2/example1/SRR30635841_1.fastq"
OUTPUT_FILE = "tests/data/sars-cov-2/example1/output.fastq"
REFERENCE_FILE = "tests/data/sars-cov-2/references/SARS-CoV-2-reference.fasta"
PRIMER_FILE = "tests/data/sars-cov-2/primers/ARTIC-V5.3.2.fasta"
AMPLICON_TYPE = "end-to-end"
COMPARE_FILE = "tests/data/sars-cov-2/example1/expected_output.fastq"


@pytest.fixture()
def setup_teardown() -> Generator[list[str], None, None]:
    args = [
        "--input",
        INPUT_FILE,
        "--output",
        OUTPUT_FILE,
        "--reference",
        REFERENCE_FILE,
        "--primer",
        PRIMER_FILE,
        "--amplicon-type",
        AMPLICON_TYPE,
    ]

    yield args

    if os.path.exists(OUTPUT_FILE):
        os.remove(OUTPUT_FILE)


def order_fastq_by_name(fastq_lines: list[str]) -> list[str]:
    chunks = [fastq_lines[i : i + 4] for i in range(0, len(fastq_lines), 4)]
    sorted_chunks = sorted(chunks, key=lambda x: x[0])
    sorted_fastq = [line for chunk in sorted_chunks for line in chunk]
    return sorted_fastq


def compare_outputs(output_file: str, expected_output_file: str) -> None:
    with open(output_file, "r") as output, open(
        expected_output_file, "r"
    ) as expected_output:
        output_lines = output.readlines()
        expected_output_lines = expected_output.readlines()
        # the fastq file is unordered, so we need to sort it by the read name
        output_lines = order_fastq_by_name(output_lines)
        expected_output_lines = order_fastq_by_name(expected_output_lines)
        for line1, line2 in zip(output_lines, expected_output_lines):
            assert line1 == line2, f"Output: {line1}, Expected Output: {line2}"


class TestAmpliGoneE2e:
    def test_ampligone(self, setup_teardown: list[str]) -> None:
        main(setup_teardown)
        compare_outputs(OUTPUT_FILE, COMPARE_FILE)
