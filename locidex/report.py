import json
import os
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from copy import deepcopy
from datetime import datetime
from dataclasses import dataclass, asdict, fields
import pandas as pd
import logging
import errno
from typing import Any
from locidex.classes.seq_intake import seq_intake
from locidex.constants import START_CODONS, STOP_CODONS, DBConfig
from locidex.utils import calc_md5
from locidex.version import __version__



logger = logging.getLogger(__name__)
logging.basicConfig(filemode=sys.stderr, level=logging.DEBUG)


@dataclass
class Parameters:
    mode: str 
    min_match_ident: str
    min_match_cov: str 
    max_ambiguous: str 
    max_internal_stops: str

@dataclass
class Data:
    sample_name: str
    profile: dict
    seq_data: dict

    def __getitem__(self, name: str) -> Any:
        return getattr(self, str(name))
    
    def __setitem__(self, key: str, value: str) -> None:
        setattr(self, key, value)

@dataclass
class ReportData:
    db_info: DBConfig
    parameters: Parameters
    data: Data

    def __getitem__(self, name: str) -> Any:
        return getattr(self, str(name))
    
    def __setitem__(self, key: str, value: str) -> None:
        setattr(self, key, value)

    @classmethod
    def fields(cls):
        return fields(cls)

    @classmethod
    def deseriealize(cls, input: dict):
        """
        Return a ReportData object from deserialized json data
        """
        return cls(db_info=DBConfig(**input["db_info"]), 
            parameters=Parameters(**input["parameters"]), 
            data=Data(**input["data"]))



def add_args(parser=None):

    if parser is None:
        parser = ArgumentParser(
            description="Locidex Report: Generate a report from search results")
    parser.add_argument('-i','--input', type=str, required=True,help='Input seq_store file to report')
    parser.add_argument('--fasta', type=str, required=False,help='Optional: Query fasta file used to generate search results')
    parser.add_argument('-c', '--config', type=str, required=False, help='Locidex parameter config file (json)')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file to put results')
    parser.add_argument('-n', '--name', type=str, required=False, help='Sample name to include default=filename')
    parser.add_argument('-m', '--mode', type=str, required=False, help='Allele profile assignment [normal,conservative,fuzzy]',default='normal')
    parser.add_argument('-p', '--prop', type=str, required=False, help='Metadata label to use for aggregation',default='locus_name')
    parser.add_argument('-a', '--max_ambig', type=int, required=False, help='Maximum number of ambiguous characters allowed in a sequence',default=0)
    parser.add_argument('-s', '--max_stop', type=int, required=False, help='Maximum number of internal stop codons allowed in a sequence',default=0)
    parser.add_argument('-r', '--match_ident', type=float, required=False, 
                        help='Report match allele if percent difference is >= this value',default=100)
    parser.add_argument('-l','--match_cov', type=float, required=False, 
                        help='Report match allele if percent coverage is >+ this value',default=100)
    parser.add_argument('--translation_table', type=int, required=False,
                        help='output directory', default=11)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    return parser



class seq_reporter:
    data_dict = {}
    profile = {}
    loci = {}
    db_seq_info = {}
    failed_seqids = set()

    def __init__(self,data_dict,method='nucleotide',mode='normal',label='locus_name',filters={},max_ambig=0,max_int_stop=0,match_ident=0):
        self.max_ambig_count = max_ambig
        self.max_int_stop_count = max_int_stop
        self.label = label
        self.method = method
        self.mode = mode
        self.match_ident = match_ident
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
            if self.mode == 'conservative':
                count_internal_stop = self.query_seq_data[seq_id]['count_internal_stop']
                start_codon = self.query_seq_data[seq_id]["start_codon"]
                stop_codon = self.query_seq_data[seq_id]["end_codon"]
                if start_codon not in START_CODONS or stop_codon not in STOP_CODONS or count_internal_stop > 0:
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
                best_hits[qid][dbtype] = []
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
                    if value >= max_bit:
                        top_ids.append(idx)
                best_hits[qid][dbtype] = top_ids

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
        self.query_best_hits = self.calc_query_best_hit()

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
            num_queries = len(loci_lookup[locus])
            if num_queries == 1 and query_hashes[0] != '-':
                assigned_loci.add(locus_name )
            elif locus_name not in loci_lookup or len(loci_lookup[locus_name]) == 0:
                assigned_loci.add(locus_name)
                self.profile[locus_name] = '-'

        loci_names_to_assign = loci_names_to_assign - assigned_loci

        profile = deepcopy(self.locus_profile)
        for locus_name in loci_names_to_assign:
            matches = loci_lookup[locus_name ]
            num_matches = len(matches)
            if num_matches <= 1:
                assigned_loci.add(locus_name)
                continue

            for qid in matches:
                if not dbtype in self.query_best_hits[qid]:
                    continue
                best_hits = self.query_best_hits[qid][dbtype]
                best_hit_names = set()
                for l in best_hits:
                    l = str(l)
                    hinfo = self.db_seq_info[l]
                    best_hit_names.add(hinfo['locus_name'])
                if len(best_hit_names) == 1:
                    continue
                if locus_name not in best_hit_names:
                    if qid in profile[locus_name]:
                        del(profile[locus_name][qid])

        self.locus_profile = profile



    def get_matching_ref_seq_info(self,qid, dbtype):
        for hit in self.query_hits[qid][dbtype]:
            hit_id = str(hit['sseqid'])
            pident = hit['pident']
            if pident < self.match_ident:
                    continue
            hinfo = self.db_seq_info[hit_id]
            hit_name = hinfo['locus_name']
            return hinfo
        return {}
    
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
                hash_value = self.query_seq_data[seq_id][key]
                if self.mode == 'fuzzy':
                    ref_seq_hitinfo = self.get_matching_ref_seq_info(seq_id, self.method)
                    if len(ref_seq_hitinfo) > 0:
                        if self.method == 'nucleotide':
                            hash_value = ref_seq_hitinfo['dna_seq_hash']  
                        elif self.method == 'protein':
                            hash_value = ref_seq_hitinfo['aa_seq_hash'] 
                            
                allele_hashes.append(hash_value)

            num_alleles = len(allele_hashes)
            if num_alleles > 1 and self.mode == 'conservative':
                allele_hashes = ['-']
            elif num_alleles > 1 and self.mode == 'normal':
                allele_hashes = calc_md5(["".join([str(x) for x in sorted(allele_hashes)])])
            elif num_alleles == 0:
                allele_hashes = ['-']
            elif self.mode == 'fuzzy':
                allele_hashes = calc_md5(["".join([str(x) for x in sorted(allele_hashes)])])
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


def run_report(config):
    
    analysis_parameters = config

    #Input Parameters
    input_file = config['input']
    outdir = config['outdir']
    label = config['prop']
    sample_name = config['name']
    force = config['force']
    mode = config['mode']
    fasta_file = config['fasta']
    max_ambig = config['max_ambig']
    max_int_stop = config['max_stop']
    match_ident = config['match_ident']
    match_cov = config['match_cov']
    translation_table = config['translation_table']


    run_data = dict()
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = analysis_parameters

    if os.path.isdir(outdir) and not force:
        logger.critical(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        raise FileExistsError(errno.EEXIST, os.strerror(errno.EEXIST), str(outdir))

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    seq_store_dict = {}
    with open(input_file ,'r') as fh:
        seq_store_dict = json.load(fh)

    if len(seq_store_dict) == 0:
        logger.critical("seq_store from file: {} is empty".format(input_file))
        raise ValueError("seq_store from file: {} is empty".format(input_file))

    if sample_name is None:
        sample_name = seq_store_dict["query_data"]["sample_name"]

    #validate the ids
    seq_data = {}
    if fasta_file is not None:
        seq_info = seq_store_dict["query_data"]["query_seq_data"]
        seq_obj = seq_intake(fasta_file, 'fasta', 'CDS', translation_table, perform_annotation=False)
        if len(seq_info) != len(seq_obj.seq_data):
            logger.critical(f'Error the supplied fasta file: {fasta_file} ({len(seq_obj.seq_data)}) seq_store file: {input_file} ({len(seq_info)}) \
                   do not have the same number of sequences. These files must be matched')
            raise ValueError(f"Supplied fasta and seq_store have different numbers of sequences: {str(fasta_file)}, {str(input_file)}")

        for i in range(0,len(seq_obj.seq_data)):
            id = str(i)
            if id not in seq_info:
                logger.critical(f'Error {id} key from fasta file not in seq_store')
                raise KeyError(f'Error {id} key from fasta file not in seq_store')
            pid_1 = seq_info[id]["seq_id"]
            pid_2 = seq_obj.seq_data[i]["seq_id"]
            if pid_1 != pid_2:
                logger.critical(f'Error seq_store key for {id}: {pid_1} mismatched to input fasta {id}: {pid_2}. These files must be matched')
                raise KeyError(f'Error seq_store key for {id}: {pid_1} mismatched to input fasta {id}: {pid_2}. These files must be matched')
            seq_data[id] = seq_obj.seq_data[i]

    allele_obj = seq_reporter(seq_store_dict, method='nucleotide', mode=mode, label=label, filters={},max_ambig=max_ambig,max_int_stop=max_int_stop,match_ident=match_ident)



    allele_obj.filter_queries()
    allele_obj.allele_assignment('nucleotide')
    allele_obj.extract_hit_data('nucleotide').to_csv(os.path.join(outdir,"nucleotide.hits.txt"),header=True,sep="\t", index=False)
    allele_obj.extract_hit_data('protein').to_csv(os.path.join(outdir, "protein.hits.txt"), header=True, sep="\t", index=False)


    profile = ReportData(
        db_info=DBConfig(**seq_store_dict["db_info"]),
        parameters= Parameters(
            mode=mode,
            min_match_ident=match_ident,
            min_match_cov=match_cov,
            max_ambiguous=max_ambig,
            max_internal_stops=max_int_stop
        ),
        data = Data(
            sample_name = sample_name,
            profile = {sample_name: allele_obj.profile},
            seq_data=seq_data
        )
    )

    
    if len(profile.data.seq_data) > 0:
        # add locus information to seq_data
        look_up = {}
        for locus_name in profile.data.profile[sample_name]:
            h = profile.data.profile[sample_name][locus_name]
            if h not in look_up:
                look_up[h] = []
            look_up[h].append(locus_name)
        
        for seq_id in profile.data.seq_data:
            h = profile.data.seq_data[seq_id]['dna_hash']
            if h in look_up:
                profile.data.seq_data[seq_id]['locus_name'] = ",".join([str(x) for x in look_up[h]])
            else:
                profile.data.seq_data[seq_id]['locus_name'] = ''


    with open(os.path.join(outdir,"report.json"),"w") as out:
        json.dump(profile,out,indent=4, default=lambda o: o.__dict__)

    run_data['result_file'] = os.path.join(outdir,"report.json")
    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    with open(os.path.join(outdir,"run.json"),'w' ) as fh:
        fh.write(json.dumps(run_data, indent=4))


def run(cmd_args=None):
    logger.info("Beginning report")
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

    run_report(config)
    logger.info("Finishing report workflow.")


# call main function
if __name__ == '__main__':
    run()



