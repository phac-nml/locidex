import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import os
from locidex.classes.fasta import parse_fasta

@pytest.fixture
def fasta_content():
    return [
        SeqRecord(Seq("ATGCGTACGTAGCTAGC"), id="gene1|123", description=""),
        SeqRecord(Seq("CGTAGCTAGCTGACGTA"), id="gene2|456", description=""),
    ]

@pytest.fixture
def fasta_file(fasta_content):
    # Create a temporary FASTA file in text mode
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as tmp:
        SeqIO.write(fasta_content, tmp, "fasta")
        tmp_path = tmp.name
    yield tmp_path
    os.unlink(tmp_path)

def test_parse_fasta_normal(fasta_file):
    parser = parse_fasta(fasta_file)
    assert len(parser.get_seqids()) == 2, "There should be two sequences"
    seq_data = parser.get_seq_by_id("gene1|123")
    assert seq_data.seq == "ATGCGTACGTAGCTAGC", "Sequence data should match the input"

def test_parse_fasta_with_nonexistent_file():
    # Assuming that the locidex errors out with File not found error if FASTA is non existant:
    with pytest.raises(FileNotFoundError):
        parse_fasta("nonexistent.fasta")

    # parser = parse_fasta("nonexistent.fasta")
    # assert not parser.status, "Parser status should be False when file does not exist"

def test_parse_fasta_with_definitions(fasta_file):
    parser = parse_fasta(fasta_file, parse_def=True, delim="|")
    seq_data = parser.get_seq_by_id("gene1|123")
    assert seq_data.gene_name == "gene1", "Gene name should be correctly parsed from definition"
    assert seq_data.seq_id == "123", "Seq ID should be correctly parsed from definition"
