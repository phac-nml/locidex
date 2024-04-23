"""
Tests for locidex format module
"""

import pytest
from locidex import format

import os
import json
import itertools
from dataclasses import dataclass

EXPECTED_DATA_OUT = "locidex/example/format_db_mlst_out"
TEST_DATA = "locidex/example/format_db_mlst_in"

@dataclass
class FormatArgs:
    input: os.PathLike
    outdir: os.PathLike
    min_len_frac: float
    max_len_frac: float
    min_ident: float
    min_match_cov: float
    translation_table: int
    not_coding: bool
    force: bool

@pytest.fixture(scope="module")
def output_directory(tmp_path_factory):
    fn = tmp_path_factory.mktemp("build")
    return str(fn)


@pytest.fixture(scope="module")
def cmd_args(output_directory):
    command = FormatArgs(
        input=TEST_DATA,
        outdir=output_directory,
        min_len_frac=0.7,
        max_len_frac=1.3,
        min_ident=80.0,
        min_match_cov=80.0,
        translation_table=11,
        not_coding=False,
        force=True
    )
    return command


def test_format(cmd_args):
    format.run(cmd_args)

def test_outputs(output_directory):
    assert sorted(os.listdir(output_directory)) == sorted(os.listdir(EXPECTED_DATA_OUT))


def test_formatted_db_content(output_directory):
    formatted_file = "locidex.txt"
    delimiter = "\t"
    actual_file = os.path.join(output_directory, formatted_file)
    expected_file = os.path.join(EXPECTED_DATA_OUT, formatted_file)
    with open(actual_file, 'r') as act, open(expected_file, 'r') as expc:
        actual_text = sorted(itertools.chain(act.readlines()).split(delimiter), key=lambda x: x[0])
        expected_text = sorted(itertools.chain(*expc.readlines()).split(delimiter), key=lambda x: x[0])
        assert hash(" ".join(actual_text)) == hash(" ".join(expected_text))

def test_format_results(output_directory):
    """
    Verify continuity of results.json fields
    """
    primary_key = "parameters"
    comp_fields = [
        "min_len_frac",
        "max_len_frac",
        "min_ident",
        "min_match_cov",
        "translation_table",
        "not_coding",
        "force"
    ]
    result_json = "results.json"
    expected = os.path.join(EXPECTED_DATA_OUT, result_json)
    actual = os.path.join(output_directory, result_json)
    with open(actual, 'r', encoding='utf8') as act, open(expected, 'r', encoding='utf8') as expc:
        act_js = json.load(act)[primary_key]
        expc_js = json.load(expc)[primary_key]
        for i in comp_fields:
            assert act_js[i] == expc_js[i]
