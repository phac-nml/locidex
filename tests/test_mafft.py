import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from locidex.classes.mafft import mafft

@pytest.fixture
def temp_fasta(tmp_path):
    # Fixture to create a temporary FASTA file with predefined sequences.
    records = [
        SeqRecord(Seq("ATGCATGCATGC"), id="seq1", description=""),
        SeqRecord(Seq("GTACGTACGTAC"), id="seq2", description="")
    ]

    # Create a temporary FASTA file
    fasta_path = tmp_path / "example.fasta"
    SeqIO.write(records, str(fasta_path), "fasta")

    # Return the path to the temporary FASTA file
    return str(fasta_path)

def test_mafft_alignment(temp_fasta):
    # Test the functionality of the mafft class using the temporary FASTA file.
    params = {"thread": "1"}
    mafft_instance = mafft(input_fasta=temp_fasta, params=params)
    assert mafft_instance.status, f"MAFFT did not execute successfully: {mafft_instance.messages}"
    alignment = mafft_instance.get_alignment()
    assert "seq1" in alignment['alignment'], "seq1 not found in alignment"
    assert "seq2" in alignment['alignment'], "seq2 not found in alignment"
    assert len(alignment['alignment']['seq1']) > 0, "Alignment for seq1 is empty"
    assert len(alignment['alignment']['seq2']) > 0, "Alignment for seq2 is empty"

def test_mafft_failure():
    # Test the mafft class with a non-existent input file to handle failure.
    params = {"thread": "1"}
    mafft_instance = mafft(input_fasta="nonexistent.fasta", params=params)
    assert not mafft_instance.status, "MAFFT did not handle non-existent file correctly"
    expected_error_substring = "Cannot open nonexistent.fasta"
    error_found = any(expected_error_substring in message for message in mafft_instance.messages)
    assert error_found, f"Expected error message substring '{expected_error_substring}' not found in messages: {mafft_instance.messages}"