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
import logging
import errno
from locidex.classes.extractor import extractor
from locidex.classes.blast import BlastSearch, FilterOptions, BlastMakeDB
from locidex.manifest import DBData
from locidex.classes.db import search_db_conf, db_config
from locidex.classes.seq_intake import seq_intake, seq_store
from locidex.constants import FILE_TYPES, BlastColumns, BlastCommands, DBConfig, DB_EXPECTED_FILES, EXTRACT_MODES, raise_file_not_found_e
from locidex.version import __version__
from locidex.classes.aligner import perform_alignment, aligner
from locidex.utils import check_db_groups, get_format

logger = logging.getLogger(__name__)
logging.basicConfig(filemode=sys.stderr, level=logging.INFO)

def add_args(parser=None):
    if parser is None:
        parser = ArgumentParser(
            description="Locidex: Extract",)

    parser.add_argument('-i','--in_fasta', type=str, required=True,help='Query assembly sequence file (fasta)')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output directory to put results')
    parser.add_argument('-n', '--name', type=str, required=False, help='Sample name to include default=filename')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--db', type=str, required=False, help='Locidex database directory')
    group.add_argument("--db_group", type=str, required=False, help="A directory of databases containing a manifest file. Requires the db_name option to be set to select the correct db")
    parser.add_argument('--db_name', type=str, required=False, help='Name of database to perform search, used when a manifest is specified as a db')
    parser.add_argument('--db_version', type=str, required=False, help='Version of database to perform search, used when a manifest is specified as a db')
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
    parser.add_argument('--mode', type=str, required=False, help='Select from the options provided',
                        default='trim', choices=EXTRACT_MODES)
    parser.add_argument('--n_threads','-t', type=int, required=False,
                        help='CPU Threads to use', default=1)
    parser.add_argument('--format', type=str, required=False,
                        help='Format of query file [genbank,fasta]')
    parser.add_argument('--translation_table', type=int, required=False,
                        help='output directory', default=11)
    parser.add_argument('-a', '--annotate', required=False, help='Perform annotation on unannotated input fasta (Do not use if you are taking in the output of extract)',
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    return parser

def write_seq_info(seq_data,out_file):
    data = {}
    for idx,entry in enumerate(seq_data):
        data[idx] = entry
    pd.DataFrame.from_dict(data,orient='index').to_csv(out_file,sep="\t",header=True,index=False)



def run_extract(config):
    # Input Parameters
    input_fasta = Path(config['in_fasta'])
    outdir = Path(config['outdir'])
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
    db_data = DBData(db_dir=db_dir)


    if not mode in EXTRACT_MODES:
        logger.critical('Provided mode for allele extraction is not valid: {}, needs to be one of ({})'.format(mode, ", ".join(EXTRACT_MODES)))
        raise ValueError('Extraction  mode is not valid: {}, needs to be one of ({})'.format(mode))

    if sample_name == None:
        sample_name = os.path.basename(input_fasta)

    config = config | db_data.config_data.to_dict() # update config data to use db data
    run_data = dict()
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = config

    format = None
    if format is None:
        format = get_format(input_fasta)
    else:
        format = format.lower()

    if format is None or format not in FILE_TYPES:
        if format is None:
            logger.critical(f'Could not guess format for {input_fasta}')
            raise ValueError("Could not guess input type for: {}".format(str(input_fasta)))
        else:
            logger.critical(f'Format for query file must be one of {list(FILE_TYPES.keys())}, you supplied {format}')
            raise ValueError(f'Format for query file must be one of {list(FILE_TYPES.keys())}, you supplied {format}')

    seq_obj = seq_intake(input_fasta, format, 'source', translation_table, perform_annotation=False,skip_trans=True)

    # Validate database is valid
    db_database_config = search_db_conf(db_dir, DB_EXPECTED_FILES, DBConfig._keys())
    if db_database_config.status == False:
        logger.critical(f'There is an issue with provided db directory: {db_dir}\n {db_database_config.messages}')
        raise ValueError("There is an issue with the provided database: {}".format(db_dir))

    metadata_path = db_database_config.meta_file_path
    metadata_obj = db_config(metadata_path, ['meta', 'info'])

    if os.path.isdir(outdir) and not force:
        logger.critical(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        raise IsADirectoryError(errno.EISDIR, os.strerror(errno.EISDIR), str(outdir))

    db_path = os.path.join(outdir, 'blast_db')

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)
    else:
        if os.path.isdir(db_path):
            shutil.rmtree(db_path)


    if not os.path.isdir(db_path):
        os.makedirs(db_path, 0o755)

    contigs_path = os.path.join(db_path,'contigs.fasta')
    seq_data = {}
    with open(contigs_path,'w') as oh:
        for idx,seq in enumerate(seq_obj.seq_data):
            seq_data[str(idx)] = {'id':str(seq.seq_id),'seq':seq.dna_seq}
            oh.write(">{}\n{}\n".format(idx,seq.dna_seq))
    del(seq_obj)
    #TODO this should probably work on more than just nucleotides
    contigs_db = BlastMakeDB(contigs_path, DBData.nucleotide_db_type(), True, contigs_path)
    contigs_db.makeblastdb()


    blast_dir_base = os.path.join(outdir, 'blast')
    if not os.path.isdir(blast_dir_base):
        os.makedirs(blast_dir_base, 0o755)

    #blast_database_paths = db_database_config.blast_paths

    blast_params = {
        'evalue': min_evalue,
        'max_target_seqs': max_target_seqs,
        'num_threads': n_threads,
        'word_size':11
    }

    nt_db = Path("{}.fasta".format(db_data.nucleotide_blast_db))
    if not nt_db.exists():
        logger.critical("Could not find file: {}".format(nt_db))
        raise_file_not_found_e(nt_db, logger)

    filter_options = {
        'evalue':  FilterOptions(min=None, max=min_evalue, include=None),
        'pident':  FilterOptions(min=min_dna_ident, max=None, include=None),
        'qcovs':   FilterOptions(min=min_dna_match_cov, max=None, include=None),
        'qcovhsp': FilterOptions(min=min_dna_match_cov, max=None, include=None),
    }

    hit_file = os.path.join(blast_dir_base, "hsps.txt")
    # TODO is this supposed to support nucleotide and amino acid?

    obj = BlastSearch(db_data=contigs_db.output_db_path, 
                    query_path=nt_db, 
                    blast_params=blast_params, 
                    blast_method=BlastCommands.blastn,
                    blast_columns=BlastColumns._fields,
                    filter_options=filter_options)
    hit_df = obj.get_blast_data(contigs_db.output_db_path, Path(hit_file))

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

    write_seq_info(exobj.seqs,os.path.join(outdir,'seq_data.txt'))
    exobj.df.to_csv(os.path.join(outdir,'filtered.hsps.txt'),header=True,sep="\t",index=False)

    nt_db_seq_obj = seq_intake(nt_db, 'fasta', 'source', translation_table, perform_annotation=False, skip_trans=True)
    nt_db_seq_data = {}
    for idx, seq in enumerate(nt_db_seq_obj.seq_data):
        nt_db_seq_data[str(seq.seq_id)] =  seq.dna_seq

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

    run_extract(config)


# call main function
if __name__ == '__main__':
    run()

