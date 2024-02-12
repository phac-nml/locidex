import json
import shutil
import pandas as pd
import os, sys, re, collections, operator, math, time,base64
from functools import partial
from mimetypes import guess_type
import gzip
from datetime import datetime
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from locidex.version import __version__
from locidex.constants import SEARCH_RUN_DATA, FILE_TYPES, BLAST_TABLE_COLS, DB_CONFIG_FIELDS,DB_EXPECTED_FILES
from locidex.classes.seq_intake import seq_intake, seq_store
from locidex.utils import write_seq_list, write_seq_dict, filter_hsps_df
from locidex.classes.blast import blast_search, parse_blast
from locidex.classes.db import search_db_conf, db_config

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


def perform_search(query_file,results_file,db_path,blast_prog,blast_params,columns):
    return blast_search(db_path,query_file,results_file,blast_params,blast_prog,columns)


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
    min_dna_ident = cmd_args.min_dna_ident
    min_aa_ident = cmd_args.min_aa_ident
    min_evalue = cmd_args.min_evalue
    min_dna_match_cov = cmd_args.min_dna_match_cov
    min_aa_match_cov = cmd_args.min_aa_match_cov
    n_threads = cmd_args.n_threads
    translation_table = cmd_args.translation_table
    force = cmd_args.force
    format = cmd_args.format
    min_dna_len = cmd_args.min_dna_len
    max_dna_len = cmd_args.max_dna_len
    min_aa_len = cmd_args.min_aa_len
    max_aa_len = cmd_args.max_aa_len
    sample_name = cmd_args.name
    perform_annotation = cmd_args.annotate
    if not perform_annotation:
        perform_annotation = False

    if sample_name == None:
        sample_name = os.path.basename(query_file)


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
    metadata_obj = db_config(metadata_path,['meta','info'])
    blast_database_paths = db_database_config.blast_paths


    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    if format is None:
        for t in FILE_TYPES:
            for ext in FILE_TYPES[t]:
                if re.search(f"{ext}$", query_file):
                    format = t
    else:
        format = format.lower()

    if format is None or format not in FILE_TYPES:
        if format is None:
            print(f'Could not guess format for {query_file}')
        else:
            print(f'Format for query file must be one of {list(FILE_TYPES.keys())}, you supplied {format}')
        sys.exit()

    seq_obj = seq_intake(query_file,format,'CDS',translation_table,perform_annotation)

    if seq_obj.status == False:
        print(f'Something went wrong parsing query file: {query_file}, please check logs and messages:\n{seq_obj.messages}')
        sys.exit()


    blast_dir_base = os.path.join(outdir,'blast')
    if not os.path.isdir(blast_dir_base):
        os.makedirs(blast_dir_base, 0o755)

    blast_params = {
        'evalue':min_evalue,
        'max_target_seqs':max_target_seqs,
        'num_threads':n_threads,
    }



    filter_options = {
        'evalue': {'min':None, 'max':min_evalue, 'include':None},
    }

    df = pd.DataFrame.from_dict(seq_obj.seq_data)
    filtered_df = df
    filtered_df['index'] = filtered_df.index.to_list()
    hit_filters = {
        'min_dna_len': min_dna_len,
        'max_dna_len': max_dna_len,
        'min_dna_ident': min_dna_ident,
        'min_dna_match_cov':min_dna_match_cov,
        'min_aa_len': min_aa_len,
        'max_aa_len': max_aa_len,
        'min_aa_ident': min_aa_ident,
        'min_aa_match_cov': min_aa_match_cov,

    }
    store_obj = seq_store(sample_name,db_database_config.config_obj.config,metadata_obj.config['meta'],seq_obj.seq_data,BLAST_TABLE_COLS,hit_filters)

    for db_label in blast_database_paths:
        label_col = 'index'
        if db_label == 'nucleotide':
            blast_prog = 'blastn'
            seq_col = 'dna_seq'
            d = os.path.join(blast_dir_base,'nucleotide')
            filter_options['pident'] = {'min':min_dna_ident, 'max':None, 'include':None}
            filter_options['qcovs'] = {'min': min_dna_match_cov, 'max': None, 'include': None}

        elif db_label == 'protein':
            blast_prog = 'blastp'
            seq_col = 'aa_seq'
            d = os.path.join(blast_dir_base, 'protein')
            filter_options['pident'] = {'min':min_aa_ident, 'max':None, 'include':None}
            filter_options['qcovs'] = {'min': min_aa_match_cov, 'max': None, 'include': None}

        if not os.path.isdir(d):
            os.makedirs(d, 0o755)
        else:
            if os.path.isfile(os.path.join(d,"queries.fasta")):
                os.remove(os.path.join(d,"queries.fasta"))
            if os.path.isfile(os.path.join(d, "hsps.txt")):
                os.remove(os.path.join(d, "hsps.txt"))

        db_path = blast_database_paths[db_label]
        create_fasta_from_df(filtered_df, label_col, seq_col, os.path.join(d,"queries.fasta"))
        perform_search(os.path.join(d,"queries.fasta"), os.path.join(d,"hsps.txt"), db_path, blast_prog, blast_params, BLAST_TABLE_COLS)
        hit_obj = parse_blast(os.path.join(d,"hsps.txt"),BLAST_TABLE_COLS,filter_options)
        hit_df = hit_obj.df
        store_obj.add_hit_data(hit_df,db_label,'qseqid')

    store_obj.filter_hits()
    store_obj.convert_profile_to_list()

    del(filtered_df)
    with open(os.path.join(outdir,"seq_store.json"),"w") as out:
        json.dump(store_obj.record,out,indent=4)


# call main function
if __name__ == '__main__':
    run()



