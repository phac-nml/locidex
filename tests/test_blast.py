import pytest, os
import pandas as pd
import locidex.classes.blast
import shutil
from locidex.classes.blast import FilterOptions
from locidex.constants import BlastColumns



PACKAGE_ROOT = os.path.dirname(locidex.__file__)

@pytest.fixture()
def blast_search_class_init(tmpdir):
    test_dir = tmpdir
    blast_search_obj = locidex.classes.blast.blast_search(input_db_path=os.path.join(PACKAGE_ROOT, "example/build_db_mlst_out/blast/nucleotide/nucleotide"), 
                        input_query_path=os.path.join(PACKAGE_ROOT, 'example/search/NC_003198.1.fasta'),
                        output_results=os.path.join(test_dir,"hsps.txt"), 
                        blast_params={'evalue': 0.0001,'max_target_seqs': 10,'num_threads': 1}, 
                        blast_method='blastn',
                        blast_columns=BlastColumns._fields)
    #blast_search_obj.run_blast()
    return test_dir, blast_search_obj


#def test_make_mlst_database(blast_search_class_init):
#    test_dir, obj = blast_search_class_init
#    #print(test_dir)
#    #blast_search_obj.input_db_path = os.path.join(PACKAGE_ROOT, "example/build_db_mlst_out/blast/nucleotide/nucleotide.fasta")
#    #blast_search_obj.input_query_path =  os.path.join(PACKAGE_ROOT, 'example/search/NC_003198.1.fasta')
#    #blast_search_obj.output_db_path = os.path.join(tmpdir,"nucleotide_mlst_database")
#    #blast_search_obj.makeblastdb(); 
#
#    assert len([file for file in os.listdir(test_dir) if "nucleotide_mlst_database" in file]) > 0
#    
#    #blast_search_obj.input_db_path = os.path.join(test_dir, "nucleotide_mlst_database") # assign new path to freshly created database to check its validity
#    #assert blast_search_obj.is_blast_db_valid() == True



def test_run_blast_on_genome_and_check_output(blast_search_class_init):
    #test_make_mlst_database(blast_search_class_init,tmpdir)
    blast_search_obj, obj = blast_search_class_init
    #blast_search_obj.input_db_path = os.path.join(tmpdir, "nucleotide_mlst_database")
    #obj.run_blast()
    #output_blast_results_path = os.path.join(tmpdir,"hsps.txt")
    output_blast_results_path = os.path.join(blast_search_obj, "hsps.txt")
    assert os.path.exists(output_blast_results_path) == True
    print(output_blast_results_path)

    with open(output_blast_results_path, "r") as fp:
        output_blast_results_file = fp.readlines()
    assert len(output_blast_results_file) == 10
    print(output_blast_results_path)
    parse_blast_obj = locidex.classes.blast.parse_blast(input_file = output_blast_results_path,
                                                        blast_columns = BlastColumns._fields,
                                                        filter_options={'bitscore': FilterOptions(min=600, max=None, include=None)})
    assert parse_blast_obj.df.shape[0] == 7
    assert parse_blast_obj.df['bitscore'].max() == 926
    assert len([item for item in parse_blast_obj.df.columns.to_list()  if item not in BlastColumns._fields]) == 0 #check if columns in df and constant are identical

