import os
from typing import Generator

import pandas as pd
import pytest

from AmpliGone.fasta2bed import main


@pytest.fixture()
def setup() -> Generator[tuple[str, str, str, str], None, None]:
    path_to_fasta = "tests/data/sars-cov-2/primers/ARTIC-V5.3.2.fasta"
    path_to_reference = "tests/data/sars-cov-2/references/SARS-CoV-2-reference.fasta"
    path_to_output = "tests/data/sars-cov-2/primers/new.bed"
    path_to_example = "tests/data/sars-cov-2/primers/SARS-CoV-2-ARTIC-V5.3.2.scheme.bed"

    yield path_to_fasta, path_to_reference, path_to_output, path_to_example

    if os.path.exists(path_to_output):
        os.remove(path_to_output)


class TestFasta2Bed:

    def compare_bed_files(self, result: str, example: str) -> None:
        res_df = pd.read_csv(result, sep="\t", header=None)
        example_df = pd.read_csv(example, sep="\t", header=None)

        # drop the names, they are not important
        res_df = res_df.drop(columns=[3])
        example_df = example_df.drop(columns=[3])

        assert res_df.equals(example_df)

    def test_fasta2bed(self, setup: tuple[str, str, str, str]) -> None:
        path_to_fasta, path_to_reference, path_to_output, path_to_example = setup
        args = [
            "--primers",
            path_to_fasta,
            "--reference",
            path_to_reference,
            "--output",
            path_to_output,
        ]
        main(args)
        assert os.path.exists(path_to_output)
        assert os.path.getsize(path_to_output) > 0
        self.compare_bed_files(path_to_output, path_to_example)
