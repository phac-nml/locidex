"""
Test manifest module
"""

import pytest
import os
import json
import shutil
from pathlib import PosixPath, Path
from locidex import manifest
from locidex.constants import DBConfig
from dataclasses import dataclass


TEST_FAIL_AUTHOR = "locidex/example/manifest_in/fails/fails_author"
TEST_FAIL_DESC = "locidex/example/manifest_in/fails/fails_name"
TEST_PASS_MULTIPLE = "locidex/example/manifest_in/passes/pass_multiple"
TEST_PASS_SINGLE = "locidex/example/manifest_in/passes/pass_single"

@dataclass
class CMDArgs:
    input: os.PathLike


def test_db_list():
    """
    List all databases in a test directory
    """
    assert manifest.check_dbs(Path(TEST_PASS_MULTIPLE)) == [PosixPath('locidex/example/manifest_in/passes/pass_multiple/pass_three_db'), 
                                                            PosixPath('locidex/example/manifest_in/passes/pass_multiple/pass_two_db'), 
                                                            PosixPath('locidex/example/manifest_in/passes/pass_multiple/pass_one_db')]

@pytest.mark.parametrize("input_dir,output",
[
    (TEST_PASS_MULTIPLE,
    [
    (PosixPath('locidex/example/manifest_in/passes/pass_multiple/pass_three_db'), 
    DBConfig(db_name='Locidex Database 3', db_version='1.0.0', db_date='04/04/2024', db_author='test1', db_desc='test1', db_num_seqs=53, is_nucl=True, is_prot=True, nucleotide_db_name='nucleotide', protein_db_name='protein')), 
    (PosixPath('locidex/example/manifest_in/passes/pass_multiple/pass_two_db'), 
    DBConfig(db_name='Locidex Database 2', db_version='1.0.0', db_date='04/04/2024', db_author='test1', db_desc='test1', db_num_seqs=53, is_nucl=True, is_prot=True, nucleotide_db_name='nucleotide', protein_db_name='protein')), 
    (PosixPath('locidex/example/manifest_in/passes/pass_multiple/pass_one_db'), 
    DBConfig(db_name='Locidex Database 1', db_version='1.0.0', db_date='04/04/2024', db_author='test1', db_desc='test1', db_num_seqs=53, is_nucl=True, is_prot=True, nucleotide_db_name='nucleotide', protein_db_name='protein'))]),
    (TEST_PASS_SINGLE,
    [(PosixPath('locidex/example/manifest_in/passes/pass_single/pass_one_db'), 
    DBConfig(db_name='Locidex Database', db_version='1.0.0', db_date='04/04/2024', db_author='test1', db_desc='test1', db_num_seqs=53, is_nucl=True, is_prot=True, nucleotide_db_name='nucleotide', protein_db_name='protein'))]),
])
def test_pass_validate_db_files(input_dir, output):
    dbs = manifest.check_dbs(Path(input_dir))
    assert manifest.validate_db_files(dbs) == output

def test_fail_validate_db_files_author(capsys):
    with pytest.raises(AttributeError):
        manifest.check_config(Path(TEST_FAIL_AUTHOR))
        assert "Config cannot have missing values: db_author is empty" == capsys.readouterr()

def test_fail_validate_db_files_description(capsys):
    with pytest.raises(AttributeError):
        manifest.check_config(Path(TEST_FAIL_AUTHOR))
        assert "Config cannot have missing values: db_desc is empty" == capsys.readouterr()


def test_create_manifest_multiple():
    output = {'Locidex Database 3': 
                {'path': 'locidex/example/manifest_in/passes/pass_multiple/pass_three_db', 
                'config': {'db_name': 'Locidex Database 3', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}, 
            'Locidex Database 2': 
                {'path': 'locidex/example/manifest_in/passes/pass_multiple/pass_two_db', 
                'config': {'db_name': 'Locidex Database 2', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}, 
            'Locidex Database 1': 
                {'path': 'locidex/example/manifest_in/passes/pass_multiple/pass_one_db', 
                'config': {'db_name': 'Locidex Database 1', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}}
    assert manifest.create_manifest(Path(TEST_PASS_MULTIPLE)) == output

def test_create_manifest_single():
    output = {'Locidex Database': 
                {'path': 'locidex/example/manifest_in/passes/pass_single/pass_one_db', 
                'config': {'db_name': 'Locidex Database', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}}
    assert manifest.create_manifest(Path(TEST_PASS_SINGLE)) == output


#def test_no_db_author(tmpdir, capsys):
#    out_dir = os.path.join(tmpdir, "build")
#    shutil.copytree(TEST_FAIL_AUTHOR, out_dir)
#    with pytest.raises(KeyError):
#        manifest.run(CMDArgs(input=out_dir))
#        assert "is missing a needed field value for db_desc, please set one" in capsys.readouterr()
#
#def test_no_db_description(tmpdir, capsys):
#    outdir = os.path.join(tmpdir, "build")
#    shutil.copytree(TEST_FAIL_DESC, outdir)
#    with pytest.raises(KeyError):
#        manifest.run(CMDArgs(input=outdir))
#        assert "is missing a needed field value for db_author, please set one" in capsys.readouterr()

