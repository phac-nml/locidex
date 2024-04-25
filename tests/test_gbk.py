import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gzip
from locidex.classes.gbk import parse_gbk

# Helper to create a GenBank file
def create_genbank_file(path, seq_records, compress=False):
    mode = 'wt'
    _open = gzip.open if compress else open
    with _open(path, mode) as f:
        for record in seq_records:
            # Ensure the record includes necessary annotations for GenBank
            record.annotations["molecule_type"] = "DNA"
        SeqIO.write(seq_records, f, "genbank")

@pytest.fixture
def genbank_file(tmp_path):
    # Create a GenBank file for testing
    file_path = tmp_path / "test.gbk"
    seqs = [SeqRecord(Seq("ATGGAATAG"), id="test_seq1", description="A simple sequence", annotations={"molecule_type": "DNA"})]
    create_genbank_file(file_path, seqs)
    return file_path

@pytest.fixture
def gzip_genbank_file(tmp_path):
    # Create a gzip-compressed GenBank file
    file_path = tmp_path / "test_compressed.gbk.gz"
    seqs = [SeqRecord(Seq("ATGGAATAG"), id="test_seq1", description="A simple sequence", annotations={"molecule_type": "DNA"})]
    create_genbank_file(file_path, seqs, compress=True)
    return file_path

def test_initialization_nonexistent_file():
    # Ensure it handles non-existent files correctly
    with pytest.raises(FileNotFoundError):
        parse_gbk("nonexistent.gbk")

def test_initialization_existing_file(genbank_file):
    # Test initialization with an existing file
    parser = parse_gbk(str(genbank_file))
    assert parser.status, "Parser should be successfully initialized with a valid file"

def test_parse_reference_gbk(genbank_file):
    # Test the parsing functionality
    parser = parse_gbk(str(genbank_file))
    parser.parse_reference_gbk()
    assert len(parser.seq_obj) > 0, "Should parse GenBank file and store data"

def test_get_acs(genbank_file):
    # Test getting accession numbers
    parser = parse_gbk(str(genbank_file))
    acs = parser.get_acs()
    assert len(acs) > 0, "Should return list of accession numbers"

def test_get_seq_by_acs(genbank_file):
    # Test fetching sequences by accession number
    parser = parse_gbk(str(genbank_file))
    seq_data = parser.get_seq_by_acs('test_seq1')
    assert seq_data, "Should return sequence data for given accession"

def test_get_feature(genbank_file):
    # Test fetching features by accession
    parser = parse_gbk(str(genbank_file))
    feature_data = parser.get_feature('test_seq1', 'source')
    assert feature_data, "Should return specific feature data for given accession"

def test_gene_prediction_gzip_fasta(gzip_genbank_file):
    # Ensure it can handle gzip compressed files
    parser = parse_gbk(str(gzip_genbank_file))
    assert parser.status, "Should be able to handle gzip compressed GenBank files"
    assert len(parser.seq_obj) > 0, "Should parse compressed GenBank file and store data"
