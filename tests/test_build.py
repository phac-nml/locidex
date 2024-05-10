""""
Tests for locidex build

"""

import pytest, warnings
from locidex import build

import os
import json
from dataclasses import dataclass

TEST_DATA = "locidex/example/build_db_mlst_in/senterica.mlst.txt"
EXPECTED_DATA_OUT = "locidex/example/build_db_mlst_out"
EXCLUDE_COMP_FILES = frozenset(["protein.pjs", "nucleotide.njs", "nucleotide.nin", "protein.pin", "config.json", "results.json"])


@dataclass
class CMDArgs:
    input_file: os.PathLike
    outdir: os.PathLike
    name: str
    db_ver: str
    db_desc: str
    author: str
    date: str
    force: bool

@pytest.fixture(scope="module")
def output_directory(tmp_path_factory):
    fn = tmp_path_factory.mktemp("build")
    return str(fn)

@pytest.fixture(scope="module")
def cmd_args(output_directory):
    command = CMDArgs(input_file=TEST_DATA, 
                    outdir=output_directory,
                    name='Locidex Database',
                    db_ver='1.0.0',
                    db_desc='test',
                    force=True,
                    author='mw',
                    date=''
                    )
    return command

def test_build(cmd_args):
    """
    Test that no errors are raised
    """
    build.run(cmd_args)
    assert len([file  for file in os.listdir(cmd_args.outdir) if "json" in file]) > 0

def test_build_outputs(output_directory):
    """
    Test initial directory structure is the same
    """
    output_data = os.listdir(output_directory)
    expected_data = os.listdir(EXPECTED_DATA_OUT)
    assert sorted(output_data) == sorted(expected_data)

def get_all_file_paths(dir):
    """
    List all files relative paths in a directory
    """
    paths = []
    for i in os.walk(dir):
        paths.extend([(i[0], g) for g in i[2]])
    return paths

@pytest.mark.parametrize("f_name,comp_fields,primary_key", [
    ("config.json", ["db_name", "db_version", "db_author", "db_desc", 
                "db_num_seqs", "is_nucl", "is_prot", "nucleotide_db_name",
                "protein_db_name"], None),
    ("results.json", ["input_file", "name", "db_ver", "db_desc", "author", "force"], "parameters")
])
def test_config_results_json(output_directory,f_name,comp_fields, primary_key):
    """Verify that config and results files outputs are the same.

    Can probably make this exclusionary fields instead of inclusion. But I want to check that fields provided still exist.
    """
    config_file_name = f_name
    actual = os.path.join(output_directory, config_file_name)
    expected = os.path.join(EXPECTED_DATA_OUT, config_file_name)
    with open(actual, 'r', encoding='utf8') as act, open(expected, 'r', encoding='utf8') as expc:
        
        act_json = json.load(act)
        expc_json = json.load(expc)

        if primary_key is not None:
            act_json = act_json[primary_key]
            expc_json = expc_json[primary_key]

        for i in comp_fields:
            assert act_json.get(i) is not None and expc_json.get(i) is not None
            assert act_json[i] == expc_json[i]


def test_validate_field_continuity(output_directory):
    """Verify that like fields between the output files are the same. The expected should match if the actual does
    otherwise this would be a time to update tests.
    """

    config_files = os.path.join(output_directory, "config.json")
    results_files = os.path.join(output_directory, "results.json")

    with open(config_files, 'r', encoding='utf8') as conf, open(results_files, 'r', encoding="utf8") as results:
        results_fields = json.load(results)["parameters"]
        conf_fields = json.load(conf)
        assert results_fields["db_ver"] == conf_fields["db_version"]
        assert results_fields["name"] == conf_fields["db_name"]
        assert results_fields["author"] == conf_fields["db_author"]
        assert results_fields["db_desc"] == conf_fields["db_desc"]
        if results_fields["date"] == '':
            warnings.warn("In results.json the date field (database) is empty")
        else:
            assert results_fields["date"] == conf_fields["db_date"]


    
@pytest.mark.skip(reason="BLAST versions provide different database structure and output results. Need to fix version to 2.15.0 and try to find better testing logic")
def test_file_outputs(output_directory):
    """
    Make sure all required blast outputs are created
    Fails as blast version needs to be set
    """
    #test_build(cmd_args)
    output_data = get_all_file_paths(output_directory)
    expected_data = get_all_file_paths(EXPECTED_DATA_OUT)
    output_files = [i[1] for i in output_data]
    expected_files = [i[1] for i in expected_data]
    assert set(output_files) == set(expected_files)

    outputs_sorted = sorted(output_data, key=lambda x: x[1])
    expected_outputs_sorted = sorted(expected_data, key=lambda x: x[1])
    assert len(outputs_sorted) == len(expected_outputs_sorted)

    for test, expected in zip(outputs_sorted, expected_outputs_sorted):
        if test[1] in EXCLUDE_COMP_FILES and expected[1] in EXCLUDE_COMP_FILES:
            continue
        test_path = os.path.join(test[0], test[1])
        expected_path = os.path.join(expected[0], expected[1])
        with open(test_path, 'rb') as t_in, open(expected_path, 'rb') as e_in:
            test_text = hash(t_in.read())
            expected_text = hash(e_in.read())
            assert test_text == expected_text