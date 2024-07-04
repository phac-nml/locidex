import pytest
import os
import locidex
#from locidex.classes.blast import blast_search, parse_blast
from locidex.classes.blast import FilterOptions, BlastSearch, BlastMakeDB
from locidex.constants import BlastColumns, BlastCommands
from locidex.classes.extractor import extractor
from locidex.classes.seq_intake import seq_intake
from locidex.classes.db import db_config
from locidex.manifest import DBData
from pathlib import Path

#could be tested via locidex extract -i ./locidex/example/search/NC_003198.1.fasta -d ./locidex/example/build_db_mlst_out/  -o tmp --force

PACKAGE_ROOT = os.path.dirname(locidex.__file__)



def blast_db_and_search(tmpdir,input_db_path):
    #blast_search_obj = blast_search(input_db_path=input_db_path, 
    #                    input_query_path=os.path.join(PACKAGE_ROOT, 'example/build_db_mlst_out/blast/nucleotide/nucleotide.fasta'),
    #                    output_results=os.path.join(tmpdir,"hsps.txt"), blast_params={'evalue': 0.0001,'max_target_seqs': 10,'num_threads': 1}, 
    #                    blast_method='blastn',
    #                    blast_columns=BlastColumns._fields,create_db=True)
    #blast_search_obj.run_blast()
    #output_blast_results_path = os.path.join(tmpdir,"hsps.txt")
    #parse_blast_obj = parse_blast(input_file = output_blast_results_path,
    #                                                    blast_columns = BlastColumns._fields,
    #                                                    filter_options={'bitscore':{'min':600, 'max':None, 'include':None}})

    blast_search_obj = BlastMakeDB(
        db_data=Path(input_db_path),
        db_type=DBData.nucleotide_db_type(),
        parse_seqids=True,
        output_db_path=Path(input_db_path))
    blast_db_path = blast_search_obj.makeblastdb()

    blast_search_obj = BlastSearch(db_data=blast_db_path, 
                query_path=Path(os.path.join(PACKAGE_ROOT, 'example/build_db_mlst_out/blast/nucleotide/nucleotide.fasta')),
                blast_params={'evalue': 0.0001,'max_target_seqs': 10,'num_threads': 1},
                blast_method=BlastCommands.blastn,
                blast_columns=BlastColumns._fields,
                filter_options={'bitscore':FilterOptions(**{'min':600, 'max':None, 'include':None})})
    
    parse_blast_obj = blast_search_obj.get_blast_data(db_path=blast_db_path, output=tmpdir.joinpath("hsps.txt"))
    
    #return parse_blast_obj

@pytest.fixture
def seq_intake_fixture():
    # Mimicking the creation of seq_data from a given input fasta file
    input_fasta = Path(PACKAGE_ROOT).joinpath('example/search/NC_003198.1.fasta')
    format = "fasta" # Adjust this based on your file type
    translation_table = 11
    seq_obj = seq_intake(input_fasta, format, 'source', translation_table, perform_annotation=False,skip_trans=True)
    return seq_obj


def test_extractor_initialization(seq_intake_fixture, tmpdir):
    db_path=os.path.join(tmpdir,"contigs.fasta")
    nt_db_test = os.path.join(PACKAGE_ROOT,'example/build_db_mlst_out/blast/nucleotide/nucleotide.fasta')
    hit_file = os.path.join(tmpdir,"hsps.txt")
    blast_params={'evalue': 0.0001, 'max_target_seqs': 10, 'num_threads': 1}
    metadata_path = os.path.join(PACKAGE_ROOT,'example/build_db_mlst_out/meta.json')
    seq_obj = seq_intake_fixture
    seq_data={}
    with open(db_path,'w') as oh:
        for idx,seq in enumerate(seq_obj.seq_data):
            seq_data[str(idx)] = {'id':str(seq.seq_id),'seq':seq.dna_seq}
            oh.write(">{}\n{}\n".format(idx,seq.dna_seq))

    blast_db = BlastMakeDB(db_path, DBData.nucleotide_db_type(), True, db_path)
    blast_db.makeblastdb()
    blast_search_obj = BlastSearch(db_data=blast_db.output_db_path,
                                query_path=nt_db_test,
                                blast_params=blast_params,
                                blast_method=BlastCommands.blastn,
                                blast_columns=BlastColumns._fields)
    
    hit_df = blast_search_obj.get_blast_data(blast_db.output_db_path, Path(hit_file))
    loci = []; metadata_obj = db_config(metadata_path, ['meta', 'info'])
    for idx,row in hit_df.iterrows():
        qid = str(row['qseqid'])
        loci.append(metadata_obj.config['meta'][qid]['locus_name'])
    hit_df['locus_name'] = loci

    extractor_instance = extractor(
        df=hit_df,
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
        extend_threshold_ratio = 0.2,
        filter_contig_breaks=True
    )
    extractor_instance.df.to_csv(os.path.join(tmpdir,'filtered.hsps.txt'),header=True,sep="\t",index=False)  

    # Verify modifications by the class's methods
    assert not extractor_instance.df.empty, "The extractor DataFrame should not be empty after initialization."

    # Assert modifications are correctly applied
    for col in ['qstart', 'qend', 'sstart', 'send']:
        assert all(extractor_instance.df[col] >= 0), f"Column {col} should have been adjusted for zero indexing."

    # Assert sequences are correctly extracted and populated in `seqs`
    assert extractor_instance.seqs, "Sequences should be extracted and the `seqs` dictionary should be populated."

    # Assert filtering based on contig boundaries if applicable
    if extractor_instance.filter_contig_breaks:
        assert not extractor_instance.df.query("`is_5prime_boundary` == False and `is_3prime_boundary`== False").empty, "Hits on contig boundaries should be filtered out when `filter_contig_breaks` is True with both 5' and 3' boundaries set to True."

    # Assert all expected columns are present after adjustments
    expected_columns = ['ext_start', 'ext_end', 'is_complete', 'is_5prime_boundary', 'is_3prime_boundary', 'reverse', 'complement', 'is_extended']
    for col in expected_columns:
        assert col in extractor_instance.df.columns, f"Expected column {col} to be present in the DataFrame."

