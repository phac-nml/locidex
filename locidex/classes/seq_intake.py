import os
import sys

from locidex.classes.gbk import parse_gbk
from locidex.classes.fasta import ParseFasta
from locidex.utils import guess_alphabet, calc_md5, six_frame_translation, slots
from locidex.classes.prodigal import gene_prediction
from locidex.constants import DNA_AMBIG_CHARS, DNA_IUPAC_CHARS, CharacterConstants, DBConfig
from typing import NamedTuple, Optional, List
from dataclasses import dataclass, asdict

@dataclass
class HitFilters:
    min_dna_len: int
    max_dna_len: int
    min_dna_ident: float
    min_dna_match_cov: float
    min_aa_len: int
    max_aa_len: int
    min_aa_ident: float
    min_aa_match_cov: float
    dna_ambig_count: int
    __slots__ = slots(__annotations__)

@dataclass
class SeqObject:
    parent_id: str
    locus_name: str
    seq_id: str
    dna_seq: str
    dna_ambig_count: int
    dna_hash: str
    dna_len: int
    aa_seq: str
    aa_hash: str
    aa_len: int
    start_codon: Optional[str] 
    end_codon: Optional[str]
    count_internal_stop: Optional[int]
    # Manually adding slots for compatibility
    __slots__ = slots(__annotations__)

    def to_dict(self) -> dict:
        return asdict(self)

class seq_intake:
    valid_types = ['genbank','gff','gtf','fasta']
    is_file_valid = ''


    def __init__(self,input_file,file_type,feat_key='CDS',translation_table=11,perform_annotation=False,num_threads=1,skip_trans=False):
        if not input_file.exists():
            raise FileNotFoundError("File {} does not exist.".format(input_file))
        
        self.input_file = input_file
        self.file_type = file_type
        self.translation_table = translation_table
        self.feat_key = feat_key
        self.skip_trans = skip_trans
        self.num_threads = num_threads
        self.prodigal_genes = []
        self.skip_trans = False
        #self.seq_data = self.process_fasta()
        self.status = True


        if file_type == 'genbank':
            self.seq_data = self.process_gbk()
            self.status = True
        elif file_type == 'fasta' and perform_annotation==True:
            self.seq_data = self.annotate_fasta(self.input_file, num_threads=self.num_threads)
        elif file_type == 'fasta' and perform_annotation==False:
            self.seq_data = self.process_fasta()
        elif file_type == 'gff': # TODO these lists do not contain all allowed file types
            self.status = False
        elif file_type == 'gtf':
            self.status = False
        else:
            raise AttributeError

        if self.status:
            self.add_codon_data()

    
    def annotate_fasta(self, input_file, num_threads):
        sobj = gene_prediction(input_file)
        sobj.predict(num_threads)
        self.prodigal_genes = sobj.genes
        seq_data = self.process_seq_hash(sobj.sequences)
        return self.process_fasta(seq_data)

    def add_codon_data(self):
        for record in self.seq_data:
            if record.aa_len == 0:
                continue
            dna_seq = record.dna_seq.lower().replace('-','')
            dna_len = len(dna_seq)
            start_codon = ''
            stop_codon = ''
            count_internal_stop = 0
            if dna_len >= 6:
                start_codon = dna_seq[0:3]
                stop_codon = dna_seq[-3:]
                count_internal_stop = record.aa_seq[:-1].count(CharacterConstants.stop_codon)
            record.start_codon = start_codon
            record.end_codon = stop_codon
            record.count_internal_stop = count_internal_stop


    def process_gbk(self) -> list[SeqObject]:
        obj = parse_gbk(self.input_file)
        seq_data = []
        if obj.status == False:
            return False
        acs = obj.get_acs()

        for a in acs:
            features = obj.get_feature(a,self.feat_key)
            for seq in features:
                s = seq['dna_seq'].lower().replace('-','')
                for char in DNA_IUPAC_CHARS:
                    s = s.replace(char,"n")

                seq_data.append(SeqObject(
                    parent_id = a,
                    locus_name = seq['gene_name'],
                    seq_id = seq['gene_name'],
                    dna_seq = s,
                    dna_ambig_count = self.count_ambig_chars(seq['dna_seq'], DNA_AMBIG_CHARS),
                    dna_hash = calc_md5([s])[0],
                    dna_len = len(seq['dna_seq']),
                    aa_seq = seq['aa_seq'],
                    aa_hash = seq['aa_hash'],
                    aa_len = len(seq['aa_seq']),
                    start_codon=None,
                    end_codon=None,
                    count_internal_stop=None,
                ))
        return seq_data

    def process_fasta(self, seq_data = []) -> list[SeqObject]:
        obj = ParseFasta(self.input_file)
        ids = obj.get_seqids()
        for id in ids:
            features = obj.get_seq_by_id(id)
            seq = features.seq
            dtype = guess_alphabet(seq)
            dna_seq = ''
            dna_hash = ''
            dna_len = 0
            aa_seq = ''
            aa_hash = ''
            aa_len = 0
            if dtype == 'dna':
                dna_seq = seq
                for char in DNA_IUPAC_CHARS:
                    dna_seq = dna_seq.replace(char,"n")
                dna_hash = calc_md5([seq])[0]
                dna_len = len(seq)
                if self.skip_trans:
                    aa_seq = ''
                else:
                    aa_seq = six_frame_translation(dna_seq,self.translation_table)[0][0]

                aa_hash = calc_md5([aa_seq])[0]
                aa_len = len(aa_seq)
            else:
                aa_seq = seq
                aa_hash = calc_md5([aa_seq])[0]
                aa_len = len(aa_seq)

            seq_data.append(SeqObject(
                parent_id=features.gene_name,
                locus_name = features.gene_name,
                seq_id = features.seq_id,
                dna_seq = dna_seq,
                dna_ambig_count = self.count_ambig_chars(dna_seq, DNA_AMBIG_CHARS),
                dna_hash = dna_hash,
                dna_len = dna_len,
                aa_seq = aa_seq,
                aa_hash = aa_hash,
                aa_len = aa_len,
                start_codon=None,
                end_codon=None,
                count_internal_stop=None,))
        return seq_data


    def process_seq_hash(self,sequences) -> list[SeqObject]:
        seq_data = []
        for id in sequences:
            seq = sequences[id]
            dtype = guess_alphabet(seq)
            dna_seq = ''
            dna_hash = ''
            dna_len = 0
            aa_seq = ''
            aa_hash = ''
            aa_len = 0
            if dtype == 'dna':
                dna_seq = seq.lower().replace('-','')
                for char in DNA_IUPAC_CHARS:
                    dna_seq = dna_seq.replace(char,"n")
                dna_hash = calc_md5([seq])[0]
                dna_len = len(seq)
                if self.skip_trans:
                    aa_seq = ''
                else:
                    aa_seq = six_frame_translation(dna_seq, self.translation_table)[0][0]
                aa_hash = calc_md5([aa_seq])[0]
                aa_len = len(aa_seq)
            else:
                aa_seq = seq.lower().replace('-','')
                aa_hash = calc_md5([aa_seq])[0]
                aa_len = len(aa_seq)
            seq_data.append(SeqObject(
                parent_id=id,
                locus_name=id,
                seq_id=id,
                dna_seq= dna_seq,
                dna_hash= dna_hash,
                dna_ambig_count=self.count_ambig_chars(dna_seq, DNA_AMBIG_CHARS),
                dna_len=dna_len,
                aa_seq=aa_seq,
                aa_hash=aa_hash,
                aa_len=aa_len,
                start_codon=None,
                end_codon=None,
                count_internal_stop=None,
            ))

        return seq_data

    def count_ambig_chars(self,seq,chars):
        count = 0
        for char in chars:
            count+= seq.count(char)
        return count


class seq_store:
    #stored_fields = ['parent_id','locus_name','seq_id','dna_hash','dna_len','aa_hash','aa_len','start_codon','stop_codon','count_internal_stop','dna_ambig_count']
    record = {
        'db_info': {},
        'db_seq_info': {},
        'query_data': {
            'sample_name':'',
            'query_seq_data': {},
            'query_hit_columns': [],
            'query_hits': {},
            "locus_profile":{}
        }
    }

    #def __init__(self,sample_name,db_config_dict,metadata_dict,query_seq_records,blast_columns,filters={},stored_fields=[]):
    def __init__(self,sample_name,db_config_dict,metadata_dict,query_seq_records,blast_columns,filters: HitFilters):
        self.sample_name = sample_name
        self.record['query_data']['sample_name'] = sample_name
        self.add_db_config(db_config_dict)
        self.add_seq_data(query_seq_records)
        self.add_db_metadata(metadata_dict)
        self.add_hit_cols(blast_columns)
        self.filters = filters
        self.record['query_data']['sample_name'] = self.sample_name


    def add_db_config(self,conf: DBConfig):
        self.record['db_info'] = conf.to_dict()

    def add_hit_cols(self,columns):
        self.record['query_hit_columns'] = columns

    def add_seq_data(self,query_seq_records: List[SeqObject]):
        for idx, v in enumerate(query_seq_records):
            self.record['query_data']['query_seq_data'][idx] = v.to_dict()

    def add_db_metadata(self,metadata_dict):
        locus_profile = {}
        for seq_index in metadata_dict:
            self.record['db_seq_info'][seq_index] = metadata_dict[seq_index]
            locus_profile[metadata_dict[seq_index]['locus_name']] = {
                    'nucleotide':set(),
                    'protein':set()
                }
        self.record['query_data']["locus_profile"] = locus_profile

    def prepopulate_hit_data(self,df,id_col):
        ids = list(df[id_col].unique())
        for i in ids:
            if str(i) not in self.record['query_data']['query_hits']:
                self.record['query_data']['query_hits'][str(i)] = {
                    'nucleotide':[],
                    'protein':[]
                }

    def add_hit_data(self,df,dbtype,id_col,prepopulate=True):
        if prepopulate:
            self.prepopulate_hit_data(df,id_col)
        for index, row in df.iterrows():
            qid = str(row[id_col])
            self.record['query_data']['query_hits'][qid][dbtype].append(row.to_dict())

    def convert_profile_to_list(self):
        for locus_name in self.record['query_data']['locus_profile']:
            for dtype in self.record['query_data']['locus_profile'][locus_name]:
                self.record['query_data']['locus_profile'][locus_name][dtype] = list(self.record['query_data']['locus_profile'][locus_name][dtype])

    def filter_hits(self):
        query_hits = self.record['query_data']['query_hits']
        for qid in query_hits:
            for dbtype in query_hits[qid]:
                filt = []
                pbitscore = 0
                pbest_hit = ''
                for hit in query_hits[qid][dbtype]:
                    hit_id = str(hit['sseqid'])

                    qlen = hit['qlen']
                    pident = hit['pident']
                    qcovs = hit['qcovs']
                    bitscore = hit['bitscore']
                    hit_name = self.record['db_seq_info'][hit_id]['locus_name']

                    hinfo = self.record["db_seq_info"][hit_id]
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
                            min_cov = self.filters['min_dna_match_cov']
                        else:
                            min_cov = hinfo["min_dna_match_cov"]
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
                            min_cov = self.filters["aa_min_cov"]
                        else:
                            min_cov = hinfo["min_aa_match_cov"]
                        if "aa_min_ident" not in hinfo:
                            min_ident = self.filters["aa_min_ident"]
                        else:
                            min_ident = hinfo["aa_min_ident"]

                    if qlen < min_len or qlen > max_len or pident < min_ident or qcovs < min_cov:
                        continue
                    self.record['query_data']["locus_profile"][hit_name][dbtype].add(qid)
                    filt.append(hit)
                    if bitscore > pbitscore:
                        pbitscore = bitscore
                        pbest_hit = hit
                query_hits[qid][dbtype] = filt






