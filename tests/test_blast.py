import pytest, os
import pandas as pd
import locidex.classes.blast
from locidex.constants import BLAST_TABLE_COLS
from locidex.classes import run_command


PACKAGE_ROOT = os.path.dirname(locidex.__file__)

@pytest.fixture
def blast_search_class_init(tmpdir):
    blast_search_obj = locidex.classes.blast.blast_search(input_db_path=None, 
                        input_query_path=os.path.join(PACKAGE_ROOT, 'example/search/NC_003198.1.fasta'),
                        output_results=os.path.join(tmpdir,"hsps.txt"), blast_params={'evalue': 0.0001,'max_target_seqs': 10,'num_threads': 1}, 
                        blast_method='blastn',
                        blast_columns=BLAST_TABLE_COLS,create_db=False)
    return blast_search_obj


def test_make_mlst_database(blast_search_class_init, tmpdir):
    blast_search_obj = blast_search_class_init
    blast_search_obj.input_db_path = os.path.join(PACKAGE_ROOT, "example/build_db_mlst_out/blast/nucleotide/nucleotide.fasta")
    blast_search_obj.input_query_path =  os.path.join(PACKAGE_ROOT, 'example/search/NC_003198.1.fasta')
    blast_search_obj.output_db_path = os.path.join(tmpdir,"nucleotide_mlst_database")
    blast_search_obj.makeblastdb(); 

    assert len([file for file in os.listdir(tmpdir) if "nucleotide_mlst_database" in file]) > 0
    
    blast_search_obj.input_db_path = os.path.join(tmpdir, "nucleotide_mlst_database") # assign new path to freshly created database to check its validity
    assert blast_search_obj.is_blast_db_valid() == True



def test_run_blast_on_genome_and_check_output(blast_search_class_init, tmpdir):
    test_make_mlst_database(blast_search_class_init,tmpdir)
    blast_search_obj = blast_search_class_init
    blast_search_obj.input_db_path = os.path.join(tmpdir, "nucleotide_mlst_database")
    blast_search_obj.run_blast()
    output_blast_results_path = os.path.join(tmpdir,"hsps.txt")
    assert os.path.exists(output_blast_results_path) == True
    with open(output_blast_results_path, "r") as fp:
        output_blast_results_file = fp.readlines()
    assert len(output_blast_results_file) == 10

    parse_blast_obj = locidex.classes.blast.parse_blast(input_file = output_blast_results_path,
                                                        blast_columns = BLAST_TABLE_COLS,
                                                        filter_options={'bitscore':{'min':600, 'max':None, 'include':None}})
    assert parse_blast_obj.df.shape[0] == 7
    assert parse_blast_obj.df['bitscore'].max() == 926
    assert len([item for item in parse_blast_obj.df.columns.to_list()  if item not in BLAST_TABLE_COLS]) == 0 #check if columns in df and constant are identical

