import os
import locidex.classes.seq_intake
from locidex.constants import BLAST_TABLE_COLS


PACKAGE_ROOT = os.path.dirname(locidex.__file__)


def seq_intake_class_init(input_file, file_type, perform_annotation):
    return locidex.classes.seq_intake.seq_intake(input_file=input_file,
                                      file_type=file_type,feat_key='CDS',translation_table=11,
                                      perform_annotation=perform_annotation,num_threads=1,skip_trans=False)

def test_read_fasta_file():
    seq_intake_object = seq_intake_class_init(input_file=os.path.join(PACKAGE_ROOT, 'example/search/NC_003198.1.fasta'),
                          file_type='fasta', perform_annotation=True)
    assert seq_intake_object.file_type == 'fasta'
    assert seq_intake_object.feat_key == 'CDS'
    assert seq_intake_object.status == True
    assert seq_intake_object.translation_table == 11
    assert seq_intake_object.valid_types == ['genbank', 'gff', 'gtf', 'fasta']
    assert len(seq_intake_object.seq_data) > 0
    assert sum([contig['aa_len'] for contig in seq_intake_object.seq_data]) > 0
    assert any([True if 'NC_003198.1' in contig['parent_id'] else False  for contig in seq_intake_object.seq_data]) == True
    assert any([True if 'NC_003198.1' in contig['locus_name'] else False  for contig in seq_intake_object.seq_data]) == True
    assert any([True if 'NC_003198.1' in contig['seq_id'] else False  for contig in seq_intake_object.seq_data]) == True
    assert len(seq_intake_object.seq_data[0]['dna_seq']) > 0
    assert seq_intake_object.seq_data[0]['dna_len'] > 0
    assert len(seq_intake_object.seq_data[0]['dna_hash']) > 0
    assert seq_intake_object.seq_data[0]['start_codon'] == 'atg'
    assert seq_intake_object.seq_data[0]['stop_codon'] == 'taa'
    assert seq_intake_object.seq_data[0]['count_internal_stop'] == 0
    assert all([ True if contig['dna_ambig_count'] == 0 else False  for contig in seq_intake_object.seq_data]) == True




    