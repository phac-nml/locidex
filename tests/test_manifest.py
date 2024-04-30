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
from typing import Dict


TEST_FAIL_AUTHOR = "locidex/example/manifest_in/fails/fails_author"
TEST_FAIL_DESC = "locidex/example/manifest_in/fails/fails_name"
TEST_PASS_MULTIPLE = "locidex/example/manifest_in/passes/pass_multiple"
TEST_PASS_SINGLE = "locidex/example/manifest_in/passes/pass_single"

@dataclass
class CMDArgs:
    input: Path


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
    (PosixPath('pass_three_db'), 
    DBConfig(db_name='Locidex Database 3', db_version='1.0.0', db_date='04/04/2024', db_author='test1', db_desc='test1', db_num_seqs=53, is_nucl=True, is_prot=True, nucleotide_db_name='nucleotide', protein_db_name='protein')), 
    (PosixPath('pass_two_db'), 
    DBConfig(db_name='Locidex Database 2', db_version='1.0.0', db_date='04/04/2024', db_author='test1', db_desc='test1', db_num_seqs=53, is_nucl=True, is_prot=True, nucleotide_db_name='nucleotide', protein_db_name='protein')), 
    (PosixPath('pass_one_db'), 
    DBConfig(db_name='Locidex Database 1', db_version='1.0.0', db_date='04/04/2024', db_author='test1', db_desc='test1', db_num_seqs=53, is_nucl=True, is_prot=True, nucleotide_db_name='nucleotide', protein_db_name='protein'))]),
    (TEST_PASS_SINGLE,
    [(PosixPath('pass_one_db'), 
    DBConfig(db_name='Locidex Database', db_version='1.0.0', db_date='04/04/2024', db_author='test1', db_desc='test1', db_num_seqs=53, is_nucl=True, is_prot=True, nucleotide_db_name='nucleotide', protein_db_name='protein'))]),
])
def test_pass_validate_db_files(input_dir, output):
    input_path = Path(input_dir)
    dbs = manifest.check_dbs(input_path)
    assert manifest.validate_db_files(dbs, input_path) == output

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
                {'path': 'pass_three_db', 
                'config': {'db_name': 'Locidex Database 3', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}, 
            'Locidex Database 2': 
                {'path': 'pass_two_db', 
                'config': {'db_name': 'Locidex Database 2', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}, 
            'Locidex Database 1': 
                {'path': 'pass_one_db', 
                'config': {'db_name': 'Locidex Database 1', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}}
    assert manifest.create_manifest(Path(TEST_PASS_MULTIPLE)) == output

def test_create_manifest_single():
    output = {'Locidex Database': 
                {'path': 'pass_one_db', 
                'config': {'db_name': 'Locidex Database', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}}
    assert manifest.create_manifest(Path(TEST_PASS_SINGLE)) == output


def test_write_manifest(tmpdir):
    outdir = tmpdir / "build"
    shutil.copytree(TEST_PASS_MULTIPLE, outdir)
    cmd_args = CMDArgs(input=outdir)
    file_out = manifest.run(cmd_args=cmd_args)
    assert file_out.exists()


def test_read_manifest(tmpdir):
    outdir = Path(tmpdir / "build")
    shutil.copytree(TEST_PASS_MULTIPLE, outdir)
    cmd_args = CMDArgs(input=outdir)
    file_out = manifest.run(cmd_args=cmd_args)
    assert file_out.exists()
    output = {'Locidex Database 3': 
                {'path': 'pass_three_db', 
                'config': {'db_name': 'Locidex Database 3', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}, 
            'Locidex Database 2': 
                {'path': 'pass_two_db', 
                'config': {'db_name': 'Locidex Database 2', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}, 
            'Locidex Database 1': 
                {'path': 'pass_one_db', 
                'config': {'db_name': 'Locidex Database 1', 'db_version': '1.0.0', 'db_date': '04/04/2024', 'db_author': 'test1', 'db_desc': 'test1', 'db_num_seqs': 53, 'is_nucl': True, 'is_prot': True, 'nucleotide_db_name': 'nucleotide', 'protein_db_name': 'protein'}}}

    manifest_data: Dict[str, manifest.ManifestItem] = manifest.read_manifest(outdir)
    for k, v in manifest_data.items():
        comp_data = output[k]
        assert v.to_dict() == comp_data
