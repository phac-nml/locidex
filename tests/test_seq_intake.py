import os, warnings
import locidex
from locidex.classes.seq_intake import SeqObject, seq_intake, seq_store, HitFilters
from locidex.constants import BlastColumns, DB_EXPECTED_FILES, DBConfig
from locidex.classes.db import search_db_conf, db_config
from collections import Counter
from dataclasses import asdict
from pathlib import Path

PACKAGE_ROOT = os.path.dirname(locidex.__file__)


def seq_intake_class_init(input_file, file_type, perform_annotation):
   #reset global class variables to avoid ambiguous results

   obj = seq_intake(input_file=Path(input_file),
                                      file_type=file_type,feat_key='CDS',translation_table=11,
                                      perform_annotation=perform_annotation,num_threads=1,skip_trans=False)
   return obj

#@pytest.mark.skip(reason="no way of currently testing this")
def test_seq_store_class():
    db_dir = os.path.join(PACKAGE_ROOT, 'example/build_db_mlst_out')
    db_database_config = search_db_conf(db_dir, DB_EXPECTED_FILES, DBConfig._keys())
    metadata_obj = db_config(db_database_config.meta_file_path, ['meta', 'info'])
    sample_name = 'NC_003198.1.fasta'

    seq_obj = seq_intake(input_file=Path(os.path.join(PACKAGE_ROOT, 'example/search/NC_003198.1.fasta')),
                          file_type='fasta', perform_annotation=False)
    
    hit_filters = HitFilters(**{'min_dna_len': 1, 'max_dna_len': 10000000, 'min_dna_ident': 80.0, 'min_dna_match_cov': 80.0, 'min_aa_len': 1, 
                   'max_aa_len': 10000000, 'min_aa_ident': 80.0, 'min_aa_match_cov': 80.0, 'dna_ambig_count': 99999999999999})
    
    seq_store_obj = seq_store(sample_name, DBConfig(**db_database_config.config_obj.config), metadata_obj.config['meta'],
                          seq_obj.seq_data, BlastColumns._fields, hit_filters)
    
    assert list(seq_store_obj.record.keys()) == ['db_info', 'db_seq_info', 'query_data', 'query_hit_columns']
    assert list(seq_store_obj.record['db_info'].keys()) == ['db_name', 'db_version', 'db_date', 'db_author', 
                                                            'db_desc', 'db_num_seqs', 'is_nucl', 'is_prot', 'nucleotide_db_name', 'protein_db_name']
    assert seq_store_obj.record['query_data']['sample_name'] == sample_name
    if len(seq_store_obj.record['query_data']['query_seq_data']) == 1:
        assert len(seq_store_obj.record['query_data']['query_seq_data']) == 1
    else:    
        warnings.warn(f"expected len(seq_store_obj.record['query_data']['query_seq_data']) == 1 but got {len(seq_store_obj.record['query_data']['query_seq_data'])}")
    
    compare_dict = seq_store_obj.record['query_data']['query_seq_data'][0]
    assert set(compare_dict.keys()) ==  set(['parent_id', 'locus_name', 'seq_id', 'dna_hash', 'dna_len', 'aa_hash', 
                                                                                    'aa_len', 'start_codon', 'end_codon', 'count_internal_stop', 'dna_ambig_count', 'dna_seq', 'aa_seq'])
    assert list(seq_store_obj.record['query_data']['locus_profile'].keys()) == ['aroC', 'dnaN', 'hemD', 'hisD', 'purE', 'sucA', 'thrA']
    assert seq_store_obj.record['query_data']['query_hit_columns'] == []
    assert seq_store_obj.record['query_data']['query_hits'] == {}


def test_read_gbk_file():
    seq_intake_object = seq_intake_class_init(input_file=os.path.join(PACKAGE_ROOT, 'example/search/NC_003198.1.gbk'),
                        file_type='genbank', perform_annotation=False)
    
    expected_orfs = 4644
    assert seq_intake_object.file_type == 'genbank'
    assert seq_intake_object.status == True
    assert seq_intake_object.translation_table == 11

    if len(seq_intake_object.seq_data) != expected_orfs:
        msg = f"Expected ORFs number is {expected_orfs} but found {len(seq_intake_object.seq_data)}! Check class init"
        warnings.warn(msg)
    else:
        assert len(seq_intake_object.seq_data) == expected_orfs

    assert seq_intake_object.seq_data[0] == SeqObject(**{'parent_id': 'NC_003198', 'locus_name': 'STY_RS00005', 'seq_id': 'STY_RS00005', 
            'dna_seq': 'atgaaccgcatcagcaccaccaccattaccaccatcaccattaccacaggtaacggtgcgggctga', 
            'dna_hash': 'f46b7aac05dba47f42391aaa5ac25edf', 'count_internal_stop': 0, 'dna_ambig_count': 0, 'start_codon': 'atg', 'end_codon': 'tga',
            'dna_len': 66, 'aa_seq': 'mnristttittitittgngag', 'aa_hash': '8b370db9e32fd0a8362c35f3535303d8', 'aa_len': 21})
    assert all([record.start_codon for record in seq_intake_object.seq_data if record.aa_len > 0]) == True
    #assert dict(Counter(['start_codon' in record if record.aa_len > 0 else False in record for record in seq_intake_object.seq_data ])) == {True: 4325, False: 319}
    assert dict(Counter([record.start_codon is not None if record.aa_len > 0 else False for record in seq_intake_object.seq_data ])) == {True: 4325, False: 319}
    

def test_read_fasta_file():
    expected_orfs=4653
    seq_intake_object = seq_intake_class_init(input_file=Path(os.path.join(PACKAGE_ROOT, 'example/search/NC_003198.1.fasta')),
                          file_type='fasta', perform_annotation=True)
    assert seq_intake_object.file_type == 'fasta'
    assert seq_intake_object.feat_key == 'CDS'
    assert seq_intake_object.status == True
    assert seq_intake_object.translation_table == 11
    assert seq_intake_object.valid_types == ['genbank', 'gff', 'gtf', 'fasta']
    if len(seq_intake_object.seq_data) != expected_orfs:
        msg = f"Expected ORFs number is {expected_orfs} but found {len(seq_intake_object.seq_data)}! Check pyrodigal and python versions."
        warnings.warn(msg)
    assert len(seq_intake_object.seq_data) > 0
    assert sum([contig.aa_len for contig in seq_intake_object.seq_data]) > 0
    assert any([True if 'NC_003198.1' in contig.parent_id else False  for contig in seq_intake_object.seq_data]) == True
    assert any([True if 'NC_003198.1' in contig.locus_name else False  for contig in seq_intake_object.seq_data]) == True
    assert any([True if 'NC_003198.1' in contig.seq_id else False  for contig in seq_intake_object.seq_data]) == True
    assert len(seq_intake_object.seq_data[0].dna_seq) > 0
    assert seq_intake_object.seq_data[0].dna_len > 0
    assert len(seq_intake_object.seq_data[0].dna_hash) > 0
    assert any([ True if contig.dna_ambig_count == 0 else False  for contig in seq_intake_object.seq_data]) == True



