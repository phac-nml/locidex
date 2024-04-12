import pytest
import os
import pandas as pd
import locidex.classes.blast
from locidex.constants import BLAST_TABLE_COLS
from locidex.classes.extractor import extractor
from locidex.classes.seq_intake import seq_intake
from locidex.classes.db import db_config, search_db_conf

PACKAGE_ROOT = os.path.dirname(locidex.__file__)


@pytest.fixture
def blast_db_and_search(tmpdir):

    ## you need to run through db, then blast, then seq intake, then you can test extractor. 
    pass

@pytest.fixture
def seq_intake_fixture(tmpdir):
    # Mimicking the creation of seq_data from a given input fasta file
    input_fasta = os.path.join(PACKAGE_ROOT, 'example/search/NC_003198.1.fasta')
    format = "fasta" # Adjust this based on your file type
    translation_table = 11
    seq_obj = seq_intake(input_fasta, format, 'CDS', translation_table, perform_annotation=False, skip_trans=True)

    seq_data = {}
    for idx, seq in enumerate(seq_obj.seq_data):
        seq_data[str(idx)] = {'id': str(seq['seq_id']), 'seq': seq['dna_seq']}
    return seq_data


def test_extractor_initialization(blast_db_and_search, seq_intake_fixture):
    mock_df = blast_db_and_search
    seq_data = seq_intake_fixture
    extractor_instance = extractor(
        df=mock_df,
        seq_data=seq_data,
        sseqid_col='sseqid',
        queryid_col='qseqid',
        qstart_col='qstart',
        qend_col='qend',
        qlen_col='qlen',
        sstart_col='sstart',
        send_col='send',
        slen_col='slen',
        sstrand_col='sstrand',
        bitscore_col='bitscore',
        overlap_thresh=1,
        extend_threshold_ratio=0.2,
        filter_contig_breaks=False
    )

    # Verify modifications by the class's methods
    assert not extractor_instance.df.empty, "The extractor DataFrame should not be empty after initialization."

    # Assert modifications are correctly applied
    for col in ['qstart', 'qend', 'sstart', 'send']:
        assert all(extractor_instance.df[col] >= 0), f"Column {col} should have been adjusted for zero indexing."

    # Assert sequences are correctly extracted and populated in `seqs`
    assert extractor_instance.seqs, "Sequences should be extracted and the `seqs` dictionary should be populated."

    # Assert filtering based on contig boundaries if applicable
    if extractor_instance.filter_contig_breaks:
        assert not extractor_instance.df[(extractor_instance.df['is_5prime_boundary']) | (extractor_instance.df['is_3prime_boundary'])].empty, "Hits on contig boundaries should be filtered out when `filter_contig_breaks` is True."

    # Assert all expected columns are present after adjustments
    expected_columns = ['ext_start', 'ext_end', 'is_complete', 'is_5prime_boundary', 'is_3prime_boundary', 'reverse', 'complement', 'is_extended']
    for col in expected_columns:
        assert col in extractor_instance.df.columns, f"Expected column {col} to be present in the DataFrame."

