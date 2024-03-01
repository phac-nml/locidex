import json
import os
import re
import shutil
import sys
from pathlib import Path
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
import numpy as np
import pandas as pd

from locidex.classes.blast import blast_search, parse_blast
from locidex.classes.db import search_db_conf, db_config
from locidex.classes.seq_intake import seq_intake, seq_store
from locidex.constants import SEARCH_RUN_DATA, FILE_TYPES, BLAST_TABLE_COLS, DB_CONFIG_FIELDS, DB_EXPECTED_FILES
from locidex.utils import write_seq_dict
from locidex.version import __version__

def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Locidex: Extract sequence features from an assembly using a locidex database",
        formatter_class=CustomFormatter)
    parser.add_argument('-i','--in_fasta', type=str, required=True,help='Query assembly sequence file (fasta)')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output directory to put results')
    parser.add_argument('-n', '--name', type=str, required=False, help='Sample name to include default=filename')
    parser.add_argument('-d', '--db', type=str, required=False, help='Locidex database directory')
    parser.add_argument('-c', '--config', type=str, required=False, help='Locidex parameter config file (json)')
    parser.add_argument('--min_evalue', type=float, required=False, help='Minumum evalue required for match',
                        default=0.0001)
    parser.add_argument('--min_dna_len', type=int, required=False, help='Global minumum query length dna',
                        default=1)
    parser.add_argument('--min_aa_len', type=int, required=False, help='Global minumum query length aa',
                        default=1)
    parser.add_argument('--max_dna_len', type=int, required=False, help='Global maximum query length dna',
                        default=10000000)
    parser.add_argument('--max_aa_len', type=int, required=False, help='Global maximum query length aa',
                        default=10000000)
    parser.add_argument('--min_dna_ident', type=float, required=False, help='Global minumum DNA percent identity required for match',
                        default=80.0)
    parser.add_argument('--min_aa_ident', type=float, required=False, help='Global minumum AA percent identity required for match',
                        default=80.0)
    parser.add_argument('--min_dna_match_cov', type=float, required=False, help='Global minumum DNA percent hit coverage identity required for match',
                        default=80.0)
    parser.add_argument('--min_aa_match_cov', type=float, required=False, help='Global minumum AA percent hit coverage identity required for match',
                        default=80.0)
    parser.add_argument('--max_target_seqs', type=int, required=False, help='Maximum number of hit seqs per query',
                        default=10)
    parser.add_argument('--n_threads','-t', type=str, required=False,
                        help='CPU Threads to use', default=1)
    parser.add_argument('--format', type=str, required=False,
                        help='Format of query file [genbank,fasta]')
    parser.add_argument('--translation_table', type=int, required=False,
                        help='output directory', default=11)
    parser.add_argument('-a', '--annotate', required=False, help='Perform annotation on unannotated input fasta',
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')

    return parser.parse_args()

class extractor:
    def __init__(self,df,seq_data,sseqid_col,queryid_col,qstart_col,qend_col,qlen_col,sstart_col,send_col,slen_col,sstrand_col):
        self.is_complete(df,qstart_col,qend_col,qlen_col)
        self.is_contig_boundary(df,sstart_col,send_col,slen_col)
        self.set_revcomp(df,sstart_col,send_col,sstrand_col)
        self.set_extraction_pos(df,sstart_col,send_col)
        self.extract_seq(df, sseqid_col, seq_data)
        pass

    def is_5prime_complete(self,df,qstart_col):
        df['is_5prime_complete'] = np.where(df[qstart_col] == 1, True, False)

    def is_3prime_complete(self,df,qend_col,qlen_col):
        df['is_3prime_complete'] = np.where(df[qend_col] == df[qlen_col] , True, False)

    def is_complete(self,df,qstart_col,qend_col,qlen_col):
        self.is_5prime_complete(df,qstart_col)
        self.is_3prime_complete(df,qend_col,qlen_col)
        df['is_complete'] = np.where(((df['is_5prime_complete'] == True) &  (df['is_3prime_complete'] == True)), True, False)

    def is_contig_boundary(self,df,sstart_col,send_col,slen_col):
        df['is_5prime_boundary'] = np.where(df[sstart_col] == 1, True, False)
        df['is_3prime_boundary'] = np.where(df[send_col] == df[slen_col], True, False)
        df['is_on_boundary'] = np.where(((df['is_5prime_boundary'] == True) & (df['is_3prime_boundary'] == True)), True, False)

    def set_revcomp(self,df,sstart_col,send_col,strand_col):
        df['reverse'] = np.where(df[sstart_col] > df[send_col], True, False)
        df['complement'] = np.where(df[strand_col] > 'minus', True, False)

    def set_extraction_pos(self,df,start_col,end_col):
        df['ext_start'] = df[start_col]
        df['ext_end'] = df[end_col]
        for idx, row in df.iterrows():
            if row['reverse'] == True:
                tmp = row['ext_end']
                end = row['ext_start']
                start = tmp
                df.loc[idx, 'ext_start'] = start
                df.loc[idx, 'ext_end'] = end

    def extract_seq(self,df,seqid_col,seq_data):
        seqs = []
        for idx,row in df.iterrows():
            seqs.append(row['ext_start'])

        #extraction needs work














def run_extract(config):
    # Input Parameters
    input_fasta = config['in_fasta']
    outdir = config['outdir']
    db_dir = config['db']
    min_dna_ident = config['min_dna_ident']
    min_evalue = config['min_evalue']
    min_dna_match_cov = config['min_dna_match_cov']
    n_threads = config['n_threads']
    translation_table = config['translation_table']
    force = config['force']
    min_dna_len = config['min_dna_len']
    sample_name = config['name']
    max_target_seqs = config['max_target_seqs']

    if sample_name == None:
        sample_name = os.path.basename(input_fasta)

    run_data = SEARCH_RUN_DATA
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = config

    format = None
    if format is None:
        for t in FILE_TYPES:
            for ext in FILE_TYPES[t]:
                if re.search(f"{ext}$", input_fasta):
                    format = t
    else:
        format = format.lower()

    if format is None or format not in FILE_TYPES:
        if format is None:
            print(f'Could not guess format for {input_fasta}')
        else:
            print(f'Format for query file must be one of {list(FILE_TYPES.keys())}, you supplied {format}')
        sys.exit()

    seq_obj = seq_intake(input_fasta, format, 'source', translation_table, perform_annotation=False,skip_trans=True)

    # Validate database is valid
    db_database_config = search_db_conf(db_dir, DB_EXPECTED_FILES, DB_CONFIG_FIELDS)
    if db_database_config.status == False:
        print(f'There is an issue with provided db directory: {db_dir}\n {db_database_config.messages}')
        sys.exit()

    metadata_path = db_database_config.meta_file_path
    metadata_obj = db_config(metadata_path, ['meta', 'info'])

    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    db_path = os.path.join(outdir, 'blast_db')
    if not os.path.isdir(db_path):
        os.makedirs(db_path, 0o755)

    db_path = os.path.join(db_path,'contigs.fasta')
    seq_data = {}
    with open(db_path,'w') as oh:
        for idx,seq in enumerate(seq_obj.seq_data):
            seq_data[idx] = {'id':seq['seq_id'],'seq':seq['dna_seq']}
            oh.write(">{}\n{}\n".format(idx,seq['dna_seq']))
    del(seq_obj)

    blast_dir_base = os.path.join(outdir, 'blast')
    if not os.path.isdir(blast_dir_base):
        os.makedirs(blast_dir_base, 0o755)

    blast_database_paths = db_database_config.blast_paths
    blast_params = {
        'evalue': min_evalue,
        'max_target_seqs': max_target_seqs,
        'num_threads': n_threads,
    }

    hit_file = os.path.join(blast_dir_base, "hsps.txt")

    obj = blast_search(input_db_path=db_path, input_query_path="{}.fasta".format(blast_database_paths['nucleotide']),
                       output_results=hit_file, blast_params=blast_params, blast_method='blastn',
                       blast_columns=BLAST_TABLE_COLS,create_db=True)

    filter_options = {
        'evalue': {'min': None, 'max': min_evalue, 'include': None},
        'pident': {'min': min_dna_ident, 'max': None, 'include': None},
        'qcovs': {'min': min_dna_match_cov, 'max': None, 'include': None}
    }

    hit_df = parse_blast(hit_file, BLAST_TABLE_COLS, filter_options).df

    extractor(hit_df,seq_data,sseqid_col='sseqid',queryid_col='qseqid',qstart_col='qstart',qend_col='qend',qlen_col='qlen',sstart_col='sstart',send_col='send',slen_col='slen',sstrand_col='sstrand')

    print(hit_df)



def run():
    cmd_args = parse_args()
    analysis_parameters = vars(cmd_args)
    config_file = cmd_args.config

    config = {}
    if config_file is not None:
        with open(config_file) as fh:
            config = json.loads(fh.read())

    for p in analysis_parameters:
        if not p in config:
            config[p] = analysis_parameters[p]

    run_extract(config)


# call main function
if __name__ == '__main__':
    run()

