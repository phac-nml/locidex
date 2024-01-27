import shutil
import pandas as pd
import os, sys, re, collections, operator, math, time,base64
from functools import partial
from mimetypes import guess_type
import gzip
import datetime
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from locidex.version import __version__
from locidex.constants import SEARCH_RUN_DATA, FILE_TYPES, BLAST_TABLE_COLS, DB_CONFIG_FIELDS,DB_EXPECTED_FILES
from locidex.classes.seq_intake import seq_intake
from locidex.utils import write_seq_list, write_seq_dict, filter_hsps_df
from locidex.classes.blast import blast_search
from locidex.classes.db import search_db_conf

def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Locidex: Advanced searching and filtering of sequence databases using query sequences",
        formatter_class=CustomFormatter)
    parser.add_argument('-q','--query', type=str, required=True,help='Query sequence file')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output directory to put results')
    parser.add_argument('-p', '--prog', type=str, required=False, help='Blast program to use [blastn, blastp, combined]',
                        default='combined')
    parser.add_argument('-d', '--db', type=str, required=False, help='Locidex database directory')
    parser.add_argument('-c', '--config', type=str, required=False, help='Locidex parameter config file (json)')
    parser.add_argument('--min_qlen', type=int, required=False, help='Minimum length for valid query sequence',default=50)
    parser.add_argument('--max_qlen', type=int, required=False, help='Maximum length for valid query sequence',default=100000)
    parser.add_argument('--max_int_stop', type=int, required=False, help='Maximum number of stop codons',
                        default=10)
    parser.add_argument('--min_evalue', type=int, required=False, help='Minumum evalue required for match',
                        default=0.0001)
    parser.add_argument('--min_hsp_len_ratio', type=float, required=False, help='Minimum length for valid query sequence',default=0.1)
    parser.add_argument('--max_hsp_len_ratio', type=float, required=False, help='Maximum length for valid query sequence',default=2)
    parser.add_argument('--n_threads', type=str, required=False,
                        help='CPU Threads to use', default=1)
    parser.add_argument('--format', type=str, required=False,
                        help='Format of query file [genbank,fasta]')
    parser.add_argument('--translation_table', type=int, required=False,
                        help='output directory', default=11)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')

    return parser.parse_args()


def perform_search(query_file,results_file,db_path,blast_prog,blast_params,columns):
    blast_search(db_path,query_file,results_file,blast_params,blast_prog,columns)

def create_fasta_from_df(df,label_col,seq_col,out_file):
    write_seq_dict(dict(zip(df[label_col].tolist(), df[seq_col])), out_file)





def run():
    cmd_args = parse_args()
    analysis_parameters = vars(cmd_args)

    #Input Parameters
    query_file = cmd_args.query
    outdir = cmd_args.outdir
    blast_program = cmd_args.prog
    db_dir = cmd_args.db
    config_file = cmd_args.config
    min_qlen = cmd_args.min_qlen
    max_qlen = cmd_args.max_qlen
    min_evalue = cmd_args.min_evalue
    min_hsp_len_ratio = cmd_args.min_hsp_len_ratio
    max_hsp_len_ratio = cmd_args.max_hsp_len_ratio
    n_threads = cmd_args.n_threads
    translation_table = cmd_args.translation_table
    force = cmd_args.force
    format = cmd_args.format
    molecule_type = cmd_args.molecule_type

    max_target_seqs = 10

    run_data = SEARCH_RUN_DATA
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)

    #Validate database is valid
    db_database_config = search_db_conf(db_dir,DB_EXPECTED_FILES,DB_CONFIG_FIELDS)
    if db_database_config.status == False:
        print(f'There is an issue with provided db directory: {db_dir}\n {db_database_config.messages}')
        sys.exit()

    metadata_path = db_database_config.meta_file_path
    blast_database_paths = db_database_config.blast_paths



    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    if format is None:
        for t in FILE_TYPES:
            for ext in FILE_TYPES[t]:
                re.search(f"{ext}$", query_file)
    else:
        format = format.lower()

    if format is None or format not in FILE_TYPES:
        if format is None:
            print(f'Could not guess format for {query_file}')
        else:
            print(f'Format for query file must be one of {list(FILE_TYPES.keys())}, you supplied {format}')
        sys.exit()

    seq_obj = seq_intake(query_file,format,'CDS',translation_table)

    if seq_obj.status == False:
        print(f'Something went wrong parsing query file: {query_file}, please check logs and messages:\n{seq_obj.messages}')
        sys.exit()


    blast_dir_base = os.path.join(outdir,'blast')
    if not os.path.isdir(blast_dir_base):
        os.makedirs(blast_dir_base, 0o755)

    blast_params = {
        'evalue':min_evalue,
        'dust':'yes',
        'max_target_seqs':max_target_seqs,
        'num_threads':n_threads,
        'outfmt':'6 {}'.format(' '.join(BLAST_TABLE_COLS))
    }

    df = pd.DataFrame.from_dict(seq_obj.seq_data)
    if molecule_type == 'dna':
        filtered_df = df[(df['dna_len'] >= min_qlen) and (df['dna_len'] <= max_qlen) ]
    else:
        filtered_df = df[(df['aa_len'] >= min_qlen) and (df['aa_len'] <= max_qlen)]

    for db_label in blast_database_paths:
        label_col = 'index'
        if db_label == 'nucleotide':
            blast_prog = 'blastn'
            seq_col = 'dna_seq'
            d = os.path.join(blast_dir_base,'nucleotide')
        elif db_label == 'protein':
            blast_prog = 'blastp'
            seq_col = 'aa_seq'
            d = os.path.join(blast_dir_base, 'protein')

        if not os.path.isdir(d):
            os.makedirs(d, 0o755)
        else:
            os.remove(os.path.join(d,"queries.fasta"))
            os.remove(os.path.join(d, "hsps.txt"))

        db_path = blast_database_paths[db_label]
        create_fasta_from_df(filtered_df, label_col, seq_col, os.path.join(d,"queries.fasta"))
        perform_search(query_file, os.path.join(d,"hsps.txt"), db_path, blast_prog, blast_params, BLAST_TABLE_COLS)

    del(filtered_df)


# call main function
if __name__ == '__main__':
    run()



