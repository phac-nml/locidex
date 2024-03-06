import json
import os
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from copy import deepcopy
from datetime import datetime

import pandas as pd

from locidex.constants import SEARCH_RUN_DATA
from locidex.utils import calc_md5
from locidex.version import __version__


def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Locidex: Advanced searching and filtering of sequence databases using query sequences",
        formatter_class=CustomFormatter)
    parser.add_argument('-i','--input', type=str, required=True,help='Input file to report')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file to put results')
    parser.add_argument('-n', '--name', type=str, required=False, help='Sample name to include default=filename')
    parser.add_argument('-m', '--mode', type=str, required=False, help='Allele profile assignment [normal,conservative]',default='normal')
    parser.add_argument('-p', '--prop', type=str, required=False, help='Metadata label to use for aggregation',default='locus_name')
    parser.add_argument('-a', '--max_ambig', type=int, required=False, help='Maximum number of ambiguous characters allowed in a sequence',default=0)
    parser.add_argument('-s', '--max_stop', type=int, required=False, help='Maximum number of internal stop codons allowed in a sequence',default=0)
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
    failed_seqids = set()

    def __init__(self,data_dict,method='nucleotide',mode='normal',label='locus_name',filters={},max_ambig=0,max_int_stop=0):
        self.max_ambig_count = max_ambig
        self.max_int_stop_count = max_int_stop
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


    def filter_queries(self):
        failed_seqids = set()
        for seq_id in self.query_seq_data:
            ambig_count = int(self.query_seq_data[seq_id]['dna_ambig_count'])
            if isinstance(self.query_seq_data[seq_id]['count_internal_stop'], int):
                stop_count = int(self.query_seq_data[seq_id]['count_internal_stop'])
            else:
                stop_count = 0

            if ambig_count > self.max_ambig_count or stop_count > self.max_int_stop_count:
                failed_seqids.add(seq_id)

        self.failed_seqids =  failed_seqids

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

    def calc_query_best_hit(self):
        best_hits = {}
        for qid in self.query_hits:
            best_hits[qid] = {}
            for dbtype in self.query_hits[qid]:
                best_hits[dbtype] = []
                hit_ids = []
                hit_bit = []
                for hit in self.query_hits[qid][dbtype]:
                    hit_ids.append(str(hit['sseqid']))
                    hit_bit.append(hit['bitscore'])
                if len(hit_bit) == 0:
                    continue
                max_bit = max(hit_bit)
                top_ids = []
                for idx, value in enumerate(hit_bit):
                    if value == max_bit:
                        top_ids.append(idx)
                best_hits[dbtype] = top_ids
        return best_hits

    def get_hit_locinames(self):
        hit_names = {}
        for qid in self.query_hits:
            hit_names[qid] = {}
            for dbtype in self.query_hits[qid]:
                hit_names[qid][dbtype] = set()
                for hit in self.query_hits[qid][dbtype]:
                    hit_id = str(hit['sseqid'])
                    hinfo = self.db_seq_info[hit_id]
                    hit_name = hinfo['locus_name']
                    hit_names[qid][dbtype].add(hit_name)
        return hit_names
    
    def get_loci_to_query_map(self,hit_names,dbtype):
        loci_lookup = {}
        for qid in hit_names:
            if not dbtype in hit_names[qid]:
                continue
            for l in hit_names[qid][dbtype]:
                if not l in loci_lookup:
                    loci_lookup[l] = []
                loci_lookup[l].append(qid)
        return loci_lookup

    def allele_assignment(self,dbtype):
        query_best_hits = self.calc_query_best_hit()
        hit_loci_names = self.get_hit_locinames()
        loci_lookup = self.get_loci_to_query_map(hit_loci_names,dbtype)

        for locus in loci_lookup:
            loci_lookup[locus] = list(set(loci_lookup[locus]) - self.failed_seqids)


        
        self.populate_profile()

        loci_names_to_assign = set(self.profile.keys())
        assigned_loci = set()

        #Fix the values of any loci where there is a single matching query or no matching queries
        for locus_name in self.profile:
            query_hashes = self.profile[locus_name].split(',')
            num_queries = len(query_hashes)
            if num_queries == 1 and query_hashes[0] != '-':
                assigned_loci.add(locus_name )
            elif locus_name not in loci_lookup or len(loci_lookup[locus_name]) == 0:
                assigned_loci.add(locus_name)

        loci_names_to_assign = loci_names_to_assign - assigned_loci

        profile = deepcopy(self.locus_profile)
        for locus_name in loci_names_to_assign:
            matches = loci_lookup[locus_name ]
            num_matches = len(matches)
            if num_matches <= 1:
                assigned_loci.add(locus_name)
                continue

            for qid in matches:
                if not dbtype in query_best_hits[qid]:
                    continue
                best_hits = query_best_hits[qid][dbtype]
                best_hit_names = set()
                for l in best_hits:
                    hinfo = self.db_seq_info[l]
                    best_hit_names.add(hinfo['locus_name'])
                if len(best_hit_names) == 1:
                    continue
                if locus_name not in best_hit_names:
                    if qid in profile[locus_name]:
                        del(profile[locus_name][qid])

        self.locus_profile = profile
        self.populate_profile()

    def populate_profile(self):
        for locus_name in self.profile:
            values = set()
            if locus_name in self.locus_profile:
                values = set(self.locus_profile[locus_name][self.method])
            allele_hashes = []
            for seq_id in values:
                if self.method == 'nucleotide':
                    key = "dna_hash"
                elif self.method == 'protein':
                    key = "aa_hash"
                allele_hashes.append(self.query_seq_data[seq_id][key])

            num_alleles = len(allele_hashes)
            if num_alleles > 1 and self.mode == 'conservative':
                allele_hashes = ['-']
            elif num_alleles > 1 and self.mode == 'normal':
                allele_hashes = calc_md5(["".join([str(x) for x in sorted(allele_hashes)])])
            elif num_alleles == 0:
                allele_hashes = ['-']
            self.profile[locus_name] = ",".join(list(set([str(x) for x in allele_hashes])))


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
    mode = cmd_args.mode
    max_ambig = cmd_args.max_ambig
    max_int_stop = cmd_args.max_stop


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

    allele_obj = seq_reporter(seq_store_dict, method='nucleotide', mode=mode, label=label, filters={},max_ambig=max_ambig,max_int_stop=max_int_stop)


    if report_format == 'profile':
        allele_obj.filter_queries()
        allele_obj.allele_assignment('nucleotide')
        profile = {sample_name: allele_obj.profile}
        with open(os.path.join(outdir,"profile.json"),"w") as out:
            json.dump(profile,out,indent=4)

    if report_format == 'profile':
        allele_obj.extract_hit_data('nucleotide').to_csv(os.path.join(outdir,"nucleotide.hits.txt"),header=True,sep="\t", index=False)
        allele_obj.extract_hit_data('protein').to_csv(os.path.join(outdir, "protein.hits.txt"), header=True, sep="\t", index=False)



# call main function
if __name__ == '__main__':
    run()



