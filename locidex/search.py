import json
import os
import re
import sys
from pathlib import Path
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
from typing import Optional
from dataclasses import dataclass
import pandas as pd
from functools import partial

#from locidex.classes.blast import blast_search, parse_blast, FilterOptions
from locidex.classes.blast2 import BlastSearch, FilterOptions
from locidex.classes.db import search_db_conf, db_config
from locidex.manifest import DBData
from locidex.classes.seq_intake import seq_intake, seq_store, HitFilters
from locidex.constants import BlastCommands, SEARCH_RUN_DATA, FILE_TYPES, BlastColumns, DB_EXPECTED_FILES, OPTION_GROUPS, DBConfig
from locidex.utils import write_seq_dict, check_db_groups, slots
from locidex.version import __version__

def add_args(parser=None):
    if parser is None:
        parser = ArgumentParser(
            description="Locidex: Advanced searching and filtering of sequence databases using query sequences",)
    parser.add_argument('-q','--query', type=str, required=True,help='Query sequence file')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output directory to put results')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--db', type=str, required=False, help='Locidex database directory')
    group.add_argument("--db_group", type=str, required=False, help="A directory of databases containing a manifest file. Requires the db_name option to be set to select the correct db")
    parser.add_argument('-n', '--name', type=str, required=False, help='Sample name to include default=filename')
    parser.add_argument('-c', '--config', type=str, required=False, help='Locidex parameter config file (json)')
    parser.add_argument('--db_name', type=str, required=False, help='Name of database to perform search, used when a manifest is specified as a db')
    parser.add_argument('--db_version', type=str, required=False, help='Version of database to perform search, used when a manifest is specified as a db')
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
    parser.add_argument('--n_threads','-t', type=int, required=False,
                        help='CPU Threads to use', default=1)
    parser.add_argument('--format', type=str, required=False,
                        help='Format of query file [genbank,fasta]')
    parser.add_argument('--translation_table', type=int, required=False,
                        help='Table to use for translation', default=11)
    parser.add_argument('-a', '--annotate', required=False, help='Perform annotation on unannotated input fasta',
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    return parser


#def perform_search(query_file,results_file,db_path,blast_prog,blast_params,columns):
#    return blast_search(db_path,query_file,results_file,blast_params,blast_prog,columns)


def create_fasta_from_df(df,label_col,seq_col,out_file):
    return write_seq_dict(dict(zip(df[label_col].tolist(), df[seq_col])), out_file)

@dataclass
class DefaultSearchOpts:
    program: str
    seq_col: str
    pident_filter: FilterOptions
    qcovs_filter: FilterOptions
    db_dir: Path
    output_dir: str


def create_outputs(output_dir: Path, db_data: DBData, blast_params: dict, configuration: DefaultSearchOpts, filtered_df: pd.DataFrame, filter_options: dict) -> pd.DataFrame:
    """
    Create outputs of blast hits
    output_dir Path: output location of search data
    configuration DefaultSearchOpts: Pararmeters passed to 'run_search' from the cli

    This function will have some needed clean up once the cli is tidied
    """
    hsps_out = "hsps.txt" #? Need to follow up on what hsps stands for
    label_col = 'index'
    query_fasta = output_dir.joinpath("queries.fasta")
    output_hsps = output_dir.joinpath(hsps_out)
    if not output_dir.exists() or not output_dir.is_dir():
        os.makedirs(output_dir, 0o755)

    output_file = create_fasta_from_df(filtered_df, label_col=label_col, seq_col=configuration.seq_col, out_file=query_fasta)
    search_data = BlastSearch(db_data, output_file, blast_params, configuration.program, BlastColumns._fields, filter_options)
    searched_df = search_data.get_blast_data(configuration.db_dir, output_hsps)
    return searched_df

def run_search(config):

    # Input Parameters
    query_file = Path(config['query'])
    outdir = Path(config['outdir'])
    db_dir = config['db']
    min_dna_ident = config['min_dna_ident']
    min_aa_ident =config['min_aa_ident']
    min_evalue = config['min_evalue']
    min_dna_match_cov = config['min_dna_match_cov']
    min_aa_match_cov = config['min_aa_match_cov']
    n_threads = config['n_threads']
    translation_table = config['translation_table']
    force = config['force']
    format = config['format']
    min_dna_len = config['min_dna_len']
    max_dna_len = config['max_dna_len']
    min_aa_len = config['min_aa_len']
    max_aa_len = config['max_aa_len']
    sample_name = config['name']
    perform_annotation = config['annotate']
    max_target_seqs = config['max_target_seqs']

    if max_count := config.get('max_ambig_count'):
        max_ambig_count = max_count
    else:
        max_ambig_count = float('inf')

    if not perform_annotation:
        perform_annotation = False

    if sample_name == None:
        sample_name = query_file.stem


    run_data = SEARCH_RUN_DATA
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = config

    db_data = DBData(db_dir=db_dir)

    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    if format is None:
        for t in FILE_TYPES:
            for ext in FILE_TYPES[t]:
                if query_file.suffix == ext:
                    format = t
    else:
        format = format.lower()

    if format is None or format not in FILE_TYPES:
        if format is None:
            print(f'Could not guess format for {query_file}')
        else:
            print(f'Format for query file must be one of {list(FILE_TYPES.keys())}, you supplied {format}')
        sys.exit()

    seq_obj = seq_intake(input_file=query_file, file_type=format, feat_key='CDS', 
                        translation_table=translation_table, perform_annotation=perform_annotation)

    if perform_annotation:
        gbk_data = []
        for idx,genes in enumerate(seq_obj.prodigal_genes):
            f = os.path.join(outdir,f"annotations_{idx}.gbk")
            with open(f,'w') as oh:
                genes.write_genbank(oh, sequence_id=f'{sample_name}_{idx}')
            gbk_data.append(Path(f).read_text())
            os.remove(f)
        f = os.path.join(outdir, f"annotations.gbk")
        with open(f, 'w') as oh:
            oh.write("\n".join([str(x) for x in gbk_data]))

    blast_params = {
        'evalue': min_evalue,
        'max_target_seqs': max_target_seqs,
        'num_threads': n_threads,
    }

    filter_options = {
        'evalue': FilterOptions(min=None, max=min_evalue, include=None)
    }

    df = pd.DataFrame.from_dict(seq_obj.seq_data)
    filtered_df = df
    filtered_df['index'] = filtered_df.index.to_list()
    hit_filters = HitFilters(
        min_dna_len = min_dna_len,
        max_dna_len=max_dna_len,
        min_dna_ident=min_dna_ident,
        min_dna_match_cov=min_dna_match_cov,
        min_aa_len=min_aa_len,
        max_aa_len=max_aa_len,
        min_aa_ident=min_aa_ident,
        min_aa_match_cov=min_aa_match_cov,
        dna_ambig_count=max_ambig_count)
    

    store_obj = seq_store(sample_name, db_data.config_data, db_data.metadata['meta'],
                        seq_obj.seq_data, BlastColumns._fields, hit_filters)

    protein_filter = DefaultSearchOpts(
        program=BlastCommands.blastp, 
        seq_col="aa_seq", 
        pident_filter=FilterOptions(min=min_aa_ident, max=None, include=None),
        qcovs_filter=FilterOptions(min=min_aa_match_cov, max=None, include=None),
        db_dir=db_data.protein_blast_db,
        output_dir=DBData.protein_name())
    
    nucleotide_filter = DefaultSearchOpts(
        program=BlastCommands.blastn, 
        seq_col="dna_seq", 
        pident_filter=FilterOptions(min=min_dna_ident, max=None, include=None),
        qcovs_filter=FilterOptions(min=min_dna_match_cov, max=None, include=None),
        db_dir=db_data.nucleotide_blast_db,
        output_dir=DBData.nucleotide_name())

    searched_hits_col = 'qseqid'
    output_creation = partial(create_outputs, output_dir=outdir, db_data=db_data, filtered_df=filtered_df, blast_params=blast_params, filter_options=filter_options)
    if db_data.nucleotide:
        searched_data = output_creation(configuration=nucleotide_filter)
        store_obj.add_hit_data(searched_data, DBData.nucleotide_name(), searched_hits_col)
    if db_data.protein:
        searched_data = output_creation(configuration=protein_filter)
        store_obj.add_hit_data(searched_data, DBData.protein_name(), searched_hits_col)

    store_obj.filter_hits()
    store_obj.convert_profile_to_list()
    run_data['result_file'] = str(outdir.joinpath("seq_store.json"))

    with open(run_data['result_file'], "w") as fh:
        fh.write(json.dumps(store_obj.record, indent=4))

    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    with open(outdir.joinpath("run.json"),'w' ) as fh:
        fh.write(json.dumps(run_data, indent=4))


def run(cmd_args=None):
    if cmd_args is None:
        parser = add_args()
        cmd_args = parser.parse_args()
    analysis_parameters = vars(cmd_args)

    analysis_parameters = check_db_groups(analysis_params=analysis_parameters, cmd_args=cmd_args)
    
    config_file = cmd_args.config

    config = {}
    if config_file is not None:
        with open(config_file) as fh:
            config = json.loads(fh.read())

    for p in analysis_parameters:
        if not p in config:
            config[p] = analysis_parameters[p]

    run_search(config)


# call main function
if __name__ == '__main__':
    run()



