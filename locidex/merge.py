import gzip
import json
import os
import re
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
from functools import partial
from mimetypes import guess_type
from multiprocessing import Pool, cpu_count
import pandas as pd
from locidex.classes.aligner import align, parse_align
from locidex.version import __version__


def add_args(parser=None):
    if parser is None:
        parser = ArgumentParser(
            description="Locidex merge: Concatonate set of input profile.json files into  a tsv table or aligned fasta")
    parser.add_argument('-i','--input', type=str, required=True,help='Input file to report', action='append', nargs='+')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file to put results')
    parser.add_argument('--n_threads','-t', type=int, required=False,
                        help='CPU Threads to use', default=1)
    parser.add_argument('--linker','-l', type=str, required=False,
                        help='Linker sequence for alignment', default='NNNNNNNNNNNNNNNNNNNN')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-a', '--align', required=False, help='Perform alignment with individual loci to produce a concatenated alignment',
                        action='store_true')
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
            encoding = guess_type(input_files[0])[1]
            _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
            with _open(input_files[0]) as f:
                for line in f:
                    file_list.append(line.rstrip())
    return file_list

def read_file_list(file_list):
    records = {}
    for f in file_list:
        if not os.path.isfile(f):
            continue
        encoding = guess_type(f)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(f) as fh:
            data = json.load(fh)
            records[data['data']['sample_name']] = data
    return records

def extract_profiles(records):
    profile = {}
    for id in records:
        for sample_name in records[id]['data']['profile']:
            profile[sample_name] = records[id]['data']['profile'][sample_name]
    return profile

def extract_seqs(records):
    seqs = {}
    for id in records:
        if not 'seq_data' in records[id]['data']:
            continue
        seqs[id] = records[id]['data']['seq_data']
    return seqs

def write_gene_fastas(seq_data,work_dir):
    d = 0
    files = {}
    for id in seq_data:
        record = seq_data[id]
        locus_name = record['locus_name']
        if 'dna_seq' in record:
            seq = record['dna_seq']
        else:
            seq = record['aa_seq']
        out_file = os.path.join(work_dir, f"{locus_name}.fas")
        if not os.path.isfile(out_file):
            oh = open(out_file,'w')
            files[locus_name] = {'file':out_file}
        else:
            oh = open(out_file,'a')
        seq_name = f'{locus_name}|{id}|{d}'
        oh.write(f'>{seq_name}\n{seq}\n')
        oh.close()
        d+=1
    return files   

def run_merge(config):
    analysis_parameters = config

    #Input Parameters
    input_files = config['input'][0]
    outdir = config['outdir']
    perform_align = config['align']
    linker_seq = config['linker']
    n_threads = config['n_threads']
    force = config['force']


    run_data = {}
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = analysis_parameters

    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    #perform merge
    file_list = get_file_list(input_files)
    records = read_file_list(file_list)

    #create profile
    df = pd.DataFrame.from_dict(extract_profiles(records), orient='index')
    df.insert(loc=0,
              column='sample_id',
              value=df.index.tolist())
    df.to_csv(os.path.join(outdir,'profile.tsv'),index=False,header=True,sep="\t")
    sample_names = list(df['sample_id'])
    del(df)
    run_data['result_file'] = os.path.join(outdir,"profile.tsv")

    #create alignment
    if perform_align and len(records) > 1:
        pass
        work_dir = os.path.join(outdir,"raw_gene_fastas")
        if not os.path.isdir(work_dir):
            os.makedirs(work_dir, 0o755)
        
        seq_data = extract_seqs(records)
        gene_files = write_gene_fastas(seq_data,work_dir)
        pool = Pool(processes=n_threads)

        results = []
        for locus_name in gene_files:
            results.append(pool.apply_async(align, args=((gene_files[locus_name]['file'],))))

        pool.close()
        pool.join()

        loci_names = list(gene_files.keys())
        alignment = {}
        for i in range(0,len(results)):
            alignment[loci_names[i]] = parse_align(results[i])
            results[i] = None
        del(results)

        out_align = os.path.join(outdir,'loci_alignment.fas')
        oh = open(out_align,'w')
        for sample_id in sample_names:
            seq = []
            for locus_name in loci_names:
                seq.append(alignment[locus_name][sample_id])
                seq.append(linker_seq)
            oh.write('>{}\n{}'.format(sample_id,"".join(seq)))
        oh.close()
        run_data['alignment_file'] = out_align
        
    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    with open(os.path.join(outdir,"run.json"),'w' ) as fh:
        fh.write(json.dumps(run_data, indent=4))




def run(cmd_args=None):
    #cmd_args = parse_args()
    if cmd_args is None:
        parser = add_args()
        cmd_args = parser.parse_args()
    analysis_parameters = vars(cmd_args)
    config_file = cmd_args.config

    config = {}
    if config_file is not None:
        with open(config_file) as fh:
            config = json.loads(fh.read())

    for p in analysis_parameters:
        if not p in config:
            config[p] = analysis_parameters[p]

    run_merge(config)


# call main function
if __name__ == '__main__':
    run()


