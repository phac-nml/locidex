import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import gzip
from locidex.classes.prodigal import gene_prediction

# Helper function to create a FASTA file
def create_fasta_file(path, seq_records, compress=False):
    mode = 'wt'
    _open = gzip.open if compress else open
    with _open(path, mode) as f:
        SeqIO.write(seq_records, f, "fasta")

@pytest.fixture
def fasta_file(tmp_path):
    # Create a FASTA file with the 16S rRNA gene from Escherichia coli
    file_path = tmp_path / "ecoli_16S_rRNA.fasta"
    ecoli_16S_seq = Seq("AGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGA")
    seqs = [SeqRecord(ecoli_16S_seq, id="Ecoli_16S_rRNA")]
    create_fasta_file(file_path, seqs)
    return file_path

@pytest.fixture
def gzip_fasta_file(tmp_path):
    file_path = tmp_path / "ecoli_16S_rRNA_compressed.fasta.gz"
    ecoli_16S_seq = Seq("AGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGA")
    seqs = [SeqRecord(ecoli_16S_seq, id="Ecoli_16S_rRNA")]
    create_fasta_file(file_path, seqs, compress=True)
    return file_path

def test_initialization_nonexistent_file():
    non_existent_filename = "nonexistent.fasta"
    gp = gene_prediction(non_existent_filename)
    assert not gp.status, "Status should be False for non-existent file"
    assert f'Error {non_existent_filename} does not exist' in gp.messages, "Proper error message should be in messages"

def test_gene_prediction_normal_fasta(fasta_file):
    predictor = gene_prediction(str(fasta_file))
    predictor.predict()
    assert predictor.status, "Should be true if gene prediction was successful"
    assert len(predictor.sequences) > 0, "Sequences should be populated"
    assert len(predictor.genes) > 0, "Sequences should be populated"

def test_gene_prediction_gzip_fasta(gzip_fasta_file):
    predictor_gzip = gene_prediction(str(gzip_fasta_file))
    predictor_gzip.predict()
    assert predictor_gzip.status, "Should be true if gene prediction was successful"
    assert len(predictor_gzip.sequences) > 0, "Sequences should be populated"
    assert len(predictor_gzip.genes) > 0, "Sequences should be populated"