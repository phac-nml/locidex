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
from locidex.classes.extractor import extractor
from locidex.classes.blast import blast_search, parse_blast
from locidex.classes.db import search_db_conf, db_config
from locidex.classes.seq_intake import seq_intake, seq_store
from locidex.constants import SEARCH_RUN_DATA, FILE_TYPES, BLAST_TABLE_COLS, DB_CONFIG_FIELDS, DB_EXPECTED_FILES, NT_SUB
from locidex.version import __version__
from locidex.classes.aligner import perform_alignment, aligner

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
    parser.add_argument('--keep_truncated', required=False, help='Keep sequences where match is broken at the end of a sequence',
                        action='store_true')
    parser.add_argument('--mode', type=str, required=False, help='(raw, trim, snps, extend)',
                        default='trim')
    parser.add_argument('--n_threads','-t', type=int, required=False,
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
    keep_truncated = config['keep_truncated']
    sample_name = config['name']
    max_target_seqs = config['max_target_seqs']
    mode = config['mode'].lower()


    if not mode in ['snps','trim','raw','extend']:
        print(f'Provided mode for allele extraction is not valid: {mode}, needs to be one of (snps, trim, extend, raw)')
        sys.exit()

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

    db_path = os.path.join(outdir, 'blast_db')

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)
    else:
        if os.path.isdir(db_path):
            shutil.rmtree(db_path)


    if not os.path.isdir(db_path):
        os.makedirs(db_path, 0o755)

    db_path = os.path.join(db_path,'contigs.fasta')
    seq_data = {}
    with open(db_path,'w') as oh:
        for idx,seq in enumerate(seq_obj.seq_data):
            seq_data[str(idx)] = {'id':str(seq['seq_id']),'seq':seq['dna_seq']}
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
    nt_db = "{}.fasta".format(blast_database_paths['nucleotide'])
    hit_file = os.path.join(blast_dir_base, "hsps.txt")

    obj = blast_search(input_db_path=db_path, input_query_path=nt_db,
                       output_results=hit_file, blast_params=blast_params, blast_method='blastn',
                       blast_columns=BLAST_TABLE_COLS,create_db=True)

    if obj.status == False:
        print("Error something went wrong, please check error messages above")
        sys.exit()

    filter_options = {
        'evalue': {'min': None, 'max': min_evalue, 'include': None},
        'pident': {'min': min_dna_ident, 'max': None, 'include': None},
        'qcovs': {'min': min_dna_match_cov, 'max': None, 'include': None}
    }

    hit_df = parse_blast(hit_file, BLAST_TABLE_COLS, filter_options).df
    hit_df['sseqid'] = hit_df['sseqid'].astype(str)
    hit_df['qseqid'] = hit_df['qseqid'].astype(str)

    #add locus data into df
    loci = []
    for idx,row in hit_df.iterrows():
        qid = str(row['qseqid'])
        loci.append(metadata_obj.config['meta'][qid]['locus_name'])

    hit_df['locus_name'] = loci

    filt_trunc = True
    if keep_truncated:
        filt_trunc = False

    exobj = extractor(hit_df,seq_data,sseqid_col='sseqid',queryid_col='qseqid',qstart_col='qstart',qend_col='qend',
                      qlen_col='qlen',sstart_col='sstart',send_col='send',slen_col='slen',sstrand_col='sstrand',
                      bitscore_col='bitscore',filter_contig_breaks=filt_trunc)

    exobj.df.to_csv(os.path.join(outdir,'filtered.hsps.txt'),header=True,sep="\t",index=False)

    nt_db_seq_obj = seq_intake(nt_db, 'fasta', 'source', translation_table, perform_annotation=False, skip_trans=True)
    nt_db_seq_data = {}
    for idx, seq in enumerate(nt_db_seq_obj.seq_data):
        nt_db_seq_data[str(seq['seq_id'])] =  seq['dna_seq']

    del(nt_db_seq_obj)

    ext_seq_data = {}
    with open(os.path.join(outdir,'raw.extracted.seqs.fasta'), 'w') as oh:
        for idx,record in enumerate(exobj.seqs):
            if min_dna_len > len(record['seq']):
                continue
            seq_id = "{}:{}:{}:{}".format(record['locus_name'],record['query_id'],record['seqid'],record['id'])
            oh.write(">{}\n{}\n".format(seq_id,record['seq']))
            ext_seq_data[seq_id] = {'locus_name':record['locus_name'],
                                'ref_id':record['query_id'],
                                'ref_seq':nt_db_seq_data[record['query_id']],
                                'ext_seq':record['seq']}
    if mode == 'trim':
        aln_obj = aligner(trim_fwd=True,trim_rev=True,ext_fwd=False, ext_rev=False,fill=False, snps_only=False)
    elif mode == 'snps':
        aln_obj = aligner(trim_fwd=True, trim_rev=True, ext_fwd=True, ext_rev=True, fill=True, snps_only=True)
    elif mode == 'extend':
        aln_obj = aligner(trim_fwd=True, trim_rev=True, ext_fwd=True, ext_rev=True, fill=False, snps_only=False)

    if mode != 'raw':
        align_dir = os.path.join(outdir, 'fastas')
        if not os.path.isdir(align_dir):
            os.makedirs(align_dir, 0o755)
        perform_alignment(ext_seq_data, align_dir, n_threads)

        with open(os.path.join(outdir, 'processed.extracted.seqs.fasta'), 'w') as oh:
            for seq_id in ext_seq_data:
                record = ext_seq_data[seq_id]
                if not 'alignment' in record:
                    continue
                ref_id = record['ref_id']
                alignment = record['alignment']
                if not seq_id in alignment or not ref_id in alignment:
                    continue
                variants = aln_obj.call_seq(seq_id,alignment[ref_id],alignment[seq_id])
                oh.write(">{}\n{}\n".format(seq_id, variants['seq'].replace('-','')))

        shutil.rmtree(align_dir)





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

