import pytest
import os
import locidex
from pathlib import Path
from locidex.manifest import DBData
from locidex.classes import blast2
from locidex.constants import BlastCommands, BlastColumns


PACKAGE_ROOT = os.path.dirname(locidex.__file__)


@pytest.fixture()
def db_data():
    db_dir = DBData(Path(PACKAGE_ROOT).joinpath("example", "build_db_mlst_out"))
    return db_dir

@pytest.fixture()
def fasta():
    return Path(PACKAGE_ROOT).joinpath("example", "search", "NC_003198.1.fasta")

def test_validate_blast_db(db_data):
    
    test_class = blast2.BlastSearch(db_data, Path("home"), dict(), BlastCommands.blastn, BlastColumns._fields, dict()) 
    test_class.validate_blast_db()


def test_blast_runs(db_data, fasta, tmpdir):
    test_class = blast2.BlastSearch(db_data, fasta, dict(), BlastCommands.blastn, BlastColumns._fields, dict())  
    out_file = "out.txt"
    stdout, stderr = test_class._run_blast(db_data.nucleotide_blast_db, tmpdir / out_file)
    with open(tmpdir / out_file, "r") as fp:
        assert len(fp.readlines()) == 30


def test_blast_runs(db_data, fasta, tmpdir):
    test_class = blast2.BlastSearch(db_data, fasta, dict(), BlastCommands.blastn, BlastColumns._fields, dict())  
    out_file = tmpdir / "out.txt"
    bd = test_class.get_blast_data(db_data.nucleotide_blast_db, out_file)
    assert len(bd) == 30