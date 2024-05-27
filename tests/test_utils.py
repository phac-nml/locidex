"""
Test the ever growing util functions
"""

import pytest
from locidex import utils
from pathlib import Path
from locidex import manifest
from argparse import Namespace
from collections import namedtuple


def test_check_db_groups_pass(monkeypatch):
    nm_group = Namespace(db_group="Db1", db_name="test_name", db_version="1.0.0")
    analysis_params = {"db_group": "Db1", "db_name": "test_name", "db_version": "1.0.0"}

    def mockreturn(*args, **kwargs):
        ret_tup = namedtuple('stuff', ["db_path"])
        ret_val = ret_tup(True)
        return ret_val
    monkeypatch.setattr(manifest, "get_manifest_db", mockreturn)
    analysis_params = utils.check_db_groups(analysis_params, nm_group)
    assert analysis_params["db"]

def test_check_db_groups_fail():
    nm_group = Namespace(db_group="Db1", db_name="test_name", db_version="1.0.0")
    analysis_params = {"db_group": "Db1", "db_name": "test_name"}

    with pytest.raises(KeyError):
        analysis_params = utils.check_db_groups(analysis_params, nm_group)


@pytest.mark.parametrize( "file_in,type",
    [
        ("test.fa", "fasta"),
        ("test.fas", "fasta"),
        ("test.ffn", "fasta"),
        ("test.fna", "fasta"),
        ("test.fasta.gz", "fasta"),
        ("test.fas.gz", "fasta"),
        ("test.fa.gz", "fasta"),
        ("test.fna.gz", "fasta"),
        ("test.ffn.gz", "fasta"),
        ("test.gbk", "genbank"),
        ("test.genbank", "genbank"),
        ("test.gbf", "genbank"),
        ("test.gbk.gz", "genbank"),
        ("test.genbank.gz", "genbank"),
        ("test.gbf.gz", "genbank"),
        ("test.gbff.gz", "genbank"),
        ("test.gbff", "genbank"),
    ])
def test_get_format(file_in, type):
    """
    test get_format function
    """
    assert utils.get_format(Path(file_in)) == type