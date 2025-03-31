import os
from typing import Generator

import pandas as pd
import pytest

from AmpliGone.fasta2bed import main


@pytest.fixture()
def setup() -> Generator[tuple[str, str, str, str], None, None]:
    path_to_fasta = "tests/data/primers/ARTIC-V5.3.2.fasta"
    path_to_reference = "tests/data/references/SARS-CoV-2-reference.fasta"
    path_to_output = "tests/data/primers/new.bed"
    path_to_example = "tests/data/primers/SARS-CoV-2-ARTIC-V5.3.2.scheme.bed"

    yield path_to_fasta, path_to_reference, path_to_output, path_to_example

    if os.path.exists(path_to_output):
        os.remove(path_to_output)


class TestFasta2Bed:
    def compare_bed_files(self, result: str, example: str) -> None:
        res_df = pd.read_csv(result, sep="\t", header=None)
        example_df = pd.read_csv(example, sep="\t", header=None)

        # drop the names [3], they are not important
        # and
        # drop the score column [4], as AmpliGone uses it to store the alignment score
        # while ARTIC (files used for testing) uses it to store the primer pool
        res_df = res_df.drop(columns=[3, 4])
        example_df = example_df.drop(columns=[3, 4])

        if not res_df.equals(example_df):
            raise AssertionError(f"{result} and {example} are not equal")

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
        if not os.path.exists(path_to_output):
            raise AssertionError(f"{path_to_output} was not created")
        if os.path.getsize(path_to_output) == 0:
            raise AssertionError(f"{path_to_output} is empty")
        self.compare_bed_files(path_to_output, path_to_example)
