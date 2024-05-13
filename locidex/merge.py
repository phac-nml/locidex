import gzip
import json
import os
import re
import sys
import errno
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
from functools import partial
from mimetypes import guess_type
from multiprocessing import Pool, cpu_count
import logging
import pandas as pd
from locidex.classes.aligner import align, parse_align
from locidex.constants import DBConfig
from locidex.report import ReportData, Data, Parameters
from locidex.version import __version__

logger = logging.getLogger(__name__)
logging.basicConfig(filemode=sys.stderr, level=logging.INFO)

def add_args(parser=None):
    """
    TODO disabling alignment until test data is prepared
    """
    if parser is None:
        parser = ArgumentParser(
            description="Locidex merge: Concatonate set of input profile.json files into  a tsv table or aligned fasta")
    parser.add_argument('-i','--input', type=str, required=True,help='Input file to report', action='append', nargs='+')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file to put results')
    #parser.add_argument('--n_threads','-t', type=int, required=False,
    #                    help='CPU Threads to use', default=1)
    #parser.add_argument('--linker','-l', type=str, required=False,
    #                    help='Linker sequence for alignment', default='NNNNNNNNNNNNNNNNNNNN')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-s', '--strict', required=False, help='Only merge data produces by the same db',
                        action='store_true')
    #parser.add_argument('-a', '--align', required=False, help='Perform alignment with individual loci to produce a concatenated alignment',
    #                    action='store_true')
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    return parser


def get_file_list(input_files):
    file_list = []
    if len(input_files) > 1:
        file_list = input_files
    else:
        if re.search(".json$", input_files[0]) or re.search(".json.gz$", input_files[0]):
            file_list = input_files
        else:
            if not os.path.isfile(input_files[0]):
                logger.critical(f'Error the supplied file {input_files[0]} does not exist')
                sys.exit(errno.ENOENT)
            encoding = guess_type(input_files[0])[1]
            _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
            with _open(input_files[0]) as f:
                for line in f:
                    line = line.rstrip()
                    if not os.path.isfile(line):
                        logger.critical(f'Error the supplied file {line} does not exist')
                        sys.exit(errno.ENOENT)
                    file_list.append(line)
    return file_list

def validate_input_file(data_in: dict, db_version: str, db_name: str, perform_validation: bool) -> tuple[ReportData, str, str]:
    """
    Validate input data for usage verifying db_versions and db_names are the same
    """

    try:
        sq_data = ReportData.deseriealize(data_in)
    except KeyError:
        logger.critical("Missing fields in configuration required fields in in reported allele file. Fields required: {}".format(ReportData.fields()))
        sys.exit()
    else:

        if db_version is not None and sq_data.db_info.db_version != db_version and perform_validation:
            logger.critical("You are attempting to merge files that were created using different database versions.")
            sys.exit()
        
        if db_name is not None and sq_data.db_info.db_name != db_name and perform_validation:
            logger.critical("You are attempting to merge files that have different names.")
            sys.exit()
    
    return sq_data, sq_data.db_info.db_version, sq_data.db_info.db_name

def check_files_exist(file_list: list[os.PathLike]) -> None:
    """
    Verify that all files to be analyzed exist
    """
    for file in file_list:
        if not os.path.isfile(file):
            logger.critical(f"Error cannot open input file {file}")
            sys.exit(errno.ENOENT)


def read_file_list(file_list,perform_validation=False):
    records = {}
    db_version = None
    db_name = None

    check_files_exist(file_list)

    for f in file_list:
        encoding = guess_type(f)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(f) as fh:
            data = json.load(fh)
            sq_data, db_version, db_name = validate_input_file(data, 
                                                            db_version=db_version, 
                                                            db_name=db_name, 
                                                            perform_validation=perform_validation)     

            sample_name = sq_data.data.sample_name
            if records.get(sq_data.data.sample_name) is None:
                records[sample_name] = sq_data
            else:
                logger.critical("Duplicate sample name detected: {}".format(sq_data.data.sample_name))
                sys.exit("Attempting to merge allele profiles with the same sample name: {}".format(sq_data.data.sample_name))
    return records

def extract_profiles(records):
    profile = {}
    for id in records:
        for sample_name in records[id].data.profile:
            if profile.get(sample_name) is not None:
                logger.critical("Sample {} already exists and will not be added.")
            profile[sample_name] = records[id].data.profile[sample_name]
    return profile

def extract_seqs(records):
    seqs = {}
    for id in records:
        seqs[id] = records[id].data.seq_data
    return seqs

def write_gene_fastas(seq_data,work_dir):
    d = 0
    files = {}
    for id in seq_data:
        record = seq_data[id]
        for seq_id in record:
            if 'locus_name' not in record[seq_id]:
                continue
            locus_name = record[seq_id]['locus_name']
            if locus_name == '':
                continue
            if 'dna_seq' in record[seq_id]:
                seq = record[seq_id]['dna_seq']
            else:
                seq = record[seq_id]['aa_seq']
            out_file = os.path.join(work_dir, f"{locus_name}.fas")
            if not os.path.isfile(out_file):
                oh = open(out_file,'w')
                files[locus_name] = {'file':out_file}
            else:
                if not locus_name in files:
                    oh = open(out_file,'w')
                    files[locus_name] = {'file':out_file}
                else:
                    oh = open(out_file,'a')
            seq_name = f'{id}'
            oh.write(f'>{seq_name}\n{seq}\n')
            oh.close()
            d+=1
    return files   

def run_merge(config):
    analysis_parameters = config

    #Input Parameters
    input_files = config['input'][0]
    outdir = config['outdir']
    ###
    # Commented out as these changes will require test data
    # perform_align = config['align']
    # linker_seq = config['linker']
    # n_threads = config['n_threads']
    ###
    force = config['force']
    validate_db = config['strict']
    if validate_db is None or validate_db == '':
        validate_db = False


    run_data = {}
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = analysis_parameters

    if os.path.isdir(outdir) and not force:
        logger.critical(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit(errno.EEXIST)

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    #perform merge
    file_list = get_file_list(input_files)
    records = read_file_list(file_list,perform_validation=validate_db)

    #create profile
    df = pd.DataFrame.from_dict(extract_profiles(records), orient='index')
    df.insert(loc=0, column='sample_id', value=df.index.tolist())
    df.to_csv(os.path.join(outdir,'profile.tsv'),index=False,header=True,sep="\t")
    
    del(df)
    run_data['result_file'] = os.path.join(outdir,"profile.tsv")
    
    ######### create alignment ###############
    # Bring this back in when test data is provided
    ############################################
    #sample_names = list(df['sample_id'])
    #loci_lengths = {}
    #invalid_loci = set()
    #if perform_align and len(records) > 1:
    #    pass
    #    work_dir = os.path.join(outdir,"raw_gene_fastas")
    #    if not os.path.isdir(work_dir):
    #        os.makedirs(work_dir, 0o755)
    #    
    #    seq_data = extract_seqs(records)
    #    gene_files = write_gene_fastas(seq_data,work_dir)
    #    del(records)
    #    del(seq_data)
    #    pool = Pool(processes=n_threads)

    #    results = []
    #    for locus_name in gene_files:
    #        results.append(pool.apply_async(align, args=((gene_files[locus_name]['file'],))))

    #    pool.close()
    #    pool.join()

    #    r = []
    #    for x in results:
    #        if isinstance(x, dict):
    #            r.append(x)
    #        else:
    #            r.append(x.get())
    #    results = r
    #    loci_names = list(gene_files.keys())
    #    alignment = {}
    #    

    #    for i in range(0,len(results)):
    #        alignment[loci_names[i]] = parse_align(results[i][0])
    #        results[i] = None
    #    del(results)

    #    
    #    for sample_id in sample_names:
    #        for locus_name in loci_names:
    #            if sample_id not in alignment[locus_name]:
    #                continue
    #            loci_lengths[locus_name] = len(alignment[locus_name][sample_id])


    #    out_align = os.path.join(outdir,'loci_alignment.fas')
    #    oh = open(out_align,'w')
    #    
    #    for sample_id in sample_names:
    #        seq = []
    #        for locus_name in loci_names:
    #            if locus_name not in loci_lengths:
    #                invalid_loci.add(locus_name)
    #                continue
    #            if sample_id in alignment[locus_name]:
    #                seq.append(alignment[locus_name][sample_id])
    #            else:
    #                seq.append(''.join(['-']*loci_lengths[locus_name]))
    #            seq.append(linker_seq)
    #        oh.write('>{}\n{}\n'.format(sample_id,"".join(seq)))
    #    oh.close()
    #    run_data['alignment_file'] = out_align

    #run_data['count_valid_loci'] = len(loci_lengths.keys())
    #run_data['count_invalid_loci'] = len(list(invalid_loci))
    #run_data['valid_loci'] = ",".join(list(loci_lengths.keys()))
    #run_data['invalid_loci'] = ",".join(list(invalid_loci))
    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    with open(os.path.join(outdir,"run.json"),'w' ) as fh:
        fh.write(json.dumps(run_data, indent=4))


def run(cmd_args=None):
    if cmd_args is None:
        parser = add_args()
        cmd_args = parser.parse_args()
    analysis_parameters = vars(cmd_args)


    config = {}
    for p in analysis_parameters:
        if not p in config:
            config[p] = analysis_parameters[p]

    run_merge(config)


# call main function
if __name__ == '__main__':
    run()


