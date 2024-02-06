import json
import pandas as pd
import os, sys, re, collections, operator, math, time,base64
from functools import partial
from mimetypes import guess_type
import gzip
from datetime import datetime
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from locidex.version import __version__
from locidex.constants import SEARCH_RUN_DATA, FILE_TYPES, BLAST_TABLE_COLS, DB_CONFIG_FIELDS,DB_EXPECTED_FILES


def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Locidex: Advanced searching and filtering of sequence databases using query sequences",
        formatter_class=CustomFormatter)
    parser.add_argument('-i','--input', type=str, required=True,help='Input file to report')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file to put results')
    parser.add_argument('-n', '--name', type=str, required=False, help='Sample name to include default=filename')
    parser.add_argument('-p', '--prop', type=str, required=False, help='Metadata label to use for aggregation',default='locus_name')
    parser.add_argument('--report_format', type=str, required=False,
                        help='Report format of parsed results [profile]',default='profile')

    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')

    return parser.parse_args()


class seq_reporter:
    data_dict = {}
    profile = {}
    loci = {}
    db_seq_info = {}
    def __init__(self,data_dict,method='nucleotide',mode='conservative',label='locus_name',filters={}):
        self.label = label
        self.method = method
        self.mode = mode
        self.data_dict = data_dict
        self.db_seq_info = self.data_dict["db_seq_info"]
        self.query_seq_data = self.data_dict["query_data"]['query_seq_data']
        self.query_hits = self.data_dict["query_data"]["query_hits"]
        self.locus_profile = self.data_dict['query_data']["locus_profile"]
        self.blast_columns = self.data_dict["query_hit_columns"]
        self.build_profile()

    def build_profile(self):
        for lid in self.data_dict["db_seq_info"]:
            locus_name = self.db_seq_info[lid]["locus_name"]
            self.loci[lid] = locus_name
            self.profile[locus_name] = []

    def filter_hits(self):
        for qid in self.query_hits:
            for dbtype in self.query_hits[qid]:
                filt = []
                for hit in self.query_hits[qid][dbtype]:
                    hit_id = str(hit['sseqid'])
                    qlen = hit['qlen']
                    pident = hit['pident']
                    qcovs = hit['qcovs']
                    bitscore = hit['bitscore']
                    hinfo = self.db_seq_info[hit_id]
                    hit_name = hinfo['locus_name']


                    if dbtype == 'nucleotide':
                        if "dna_min_len" not in hinfo:
                            min_len = self.filters["dna_min_len"]
                        else:
                            min_len = hinfo["dna_min_len"]
                        if "dna_max_len" not in hinfo:
                            max_len = self.filters["dna_min_len"]
                        else:
                            max_len = hinfo["dna_max_len"]
                        if "min_dna_match_cov" not in hinfo:
                            min_cov = self.filters["min_dna_match_cov"]
                        else:
                            min_cov = hinfo["dna_min_cov"]
                        if "dna_min_ident" not in hinfo:
                            min_ident = self.filters["dna_min_ident"]
                        else:
                            min_ident = hinfo["dna_min_ident"]
                    else:
                        if qlen < hinfo["aa_min_len"] or qlen > hinfo["aa_max_len"] or pident < hinfo["aa_min_ident"]:
                            continue
                        if "aa_min_len" not in hinfo:
                            min_len = self.filters["aa_min_len"]
                        else:
                            min_len = hinfo["aa_min_len"]
                        if "aa_max_len" not in hinfo:
                            max_len = self.filters["aa_min_len"]
                        else:
                            max_len = hinfo["aa_max_len"]
                        if "min_aa_match_cov" not in hinfo:
                            min_cov = self.filters["min_aa_match_cov"]
                        else:
                            min_cov = hinfo["aa_min_cov"]
                        if "aa_min_ident" not in hinfo:
                            min_ident = self.filters["aa_min_ident"]
                        else:
                            min_ident = hinfo["aa_min_ident"]

                    if qlen < min_len or qlen > max_len or pident < min_ident or qcovs < min_cov:
                        continue
                    self.record['query_data']["locus_profile"][hit_name][dbtype].append(qid)
                    filt.append(hit)

                self.query_hits[qid][dbtype] = filt

    def populate_profile(self):
        for loci_name in self.profile:
            values = set()
            if loci_name in self.locus_profile:
                values = set(self.locus_profile[loci_name][self.method])
            allele_hashes = []
            for seq_id in values:
                if self.method == 'nucleotide':
                    key = "dna_hash"
                elif self.method == 'protein':
                    key = "aa_hash"
                allele_hashes.append(self.query_seq_data[seq_id][key])
            num_alleles = len(allele_hashes)
            if num_alleles > 1 and self.mode == 'conservative':
                allele_hashes = []
            self.profile[loci_name] = ",".join([str(x) for x in allele_hashes])


    def extract_hit_data(self,dbtype):
        data = []
        for qid in self.query_hits:
            query_name = self.query_seq_data[qid][self.label]
            for hit in self.query_hits[qid][dbtype]:
                hit_id = str(hit['sseqid'])
                hinfo = self.db_seq_info[hit_id]
                record = {'query_name':query_name}
                for fieldname in self.blast_columns:
                    record[fieldname] = hit[fieldname]
                for fieldname in hinfo:
                    record[fieldname] = hinfo[fieldname]

                data.append(record)
        return pd.DataFrame.from_dict(data)


def run():
    cmd_args = parse_args()
    analysis_parameters = vars(cmd_args)

    #Input Parameters
    input_file = cmd_args.input
    outdir = cmd_args.outdir
    label = cmd_args.prop
    report_format = cmd_args.report_format
    sample_name = cmd_args.name
    force = cmd_args.force


    if sample_name is None:
        sample_name = '.'.join(os.path.basename(input_file).split('.')[:-1])

    run_data = SEARCH_RUN_DATA
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)

    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    seq_store_dict = {}
    with open(input_file ,'r') as fh:
        seq_store_dict = json.load(fh)

    if len(seq_store_dict) == 0:
        sys.exit()

    allele_obj = seq_reporter(seq_store_dict, method='nucleotide', mode='conservative', label=label, filters={})
    if report_format == 'profile':
        allele_obj.populate_profile()
        profile = {sample_name: allele_obj.profile}
        with open(os.path.join(outdir,"profile.json"),"w") as out:
            json.dump(profile,out,indent=4)

    if report_format == 'profile':
        allele_obj.extract_hit_data('nucleotide').to_csv(os.path.join(outdir,"nucleotide.hits.txt"),header=True,sep="\t", index=False)
        allele_obj.extract_hit_data('protein').to_csv(os.path.join(outdir, "protein.hits.txt"), header=True, sep="\t", index=False)



# call main function
if __name__ == '__main__':
    run()



