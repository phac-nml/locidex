import os
import sys

from locidex.classes.gbk import parse_gbk
from locidex.classes.fasta import parse_fasta
from locidex.utils import guess_alphabet, calc_md5, six_frame_translation
from locidex.classes.prodigal import gene_prediction
from locidex.constants import DNA_AMBIG_CHARS, DNA_IUPAC_CHARS

class seq_intake:
    input_file = ''
    valid_types = ['genbank','gff','gtf','fasta']
    file_type = None
    feat_key = 'CDS'
    translation_table = 11
    is_file_valid = ''
    status = True
    messages = []
    seq_data = []
    prodigal_genes = []
    skip_trans = False

    def __init__(self,input_file,file_type,feat_key='CDS',translation_table=11,perform_annotation=False,num_threads=1,skip_trans=False):
        self.input_file = input_file
        self.file_type = file_type
        self.translation_table = translation_table
        self.feat_key = feat_key
        self.skip_trans = skip_trans

        if not os.path.isfile(self.input_file):
            self.messages.append(f'Error {self.input_file} does not exist')
            self.status = False

        if not self.status:
            return

        if file_type == 'genbank':
            self.status = self.process_gbk()
        elif file_type == 'fasta' and perform_annotation==True:
            sobj = gene_prediction(self.input_file)
            sobj.predict(num_threads)
            self.prodigal_genes = sobj.genes
            self.process_seq_hash(sobj.sequences)
            self.process_fasta()
        elif file_type == 'fasta' and perform_annotation==False:
            self.process_fasta()
        elif file_type == 'gff':
            self.status = False
        elif file_type == 'gtf':
            self.status = False
        if self.status:
            self.add_codon_data()


    def add_codon_data(self):
        for record in self.seq_data:
            if record['aa_len'] == 0:
                continue
            dna_seq = record['dna_seq'].lower().replace('-','')
            dna_len = len(dna_seq)
            start_codon = ''
            stop_codon = ''
            count_internal_stop = 0
            if dna_len >= 6:
                start_codon = dna_seq[0:3]
                stop_codon = dna_seq[-3:]
                count_internal_stop = record['aa_seq'][:-1].count('*')
            record['start_codon'] = start_codon
            record['stop_codon'] = stop_codon
            record['count_internal_stop'] = count_internal_stop


    def process_gbk(self):
        obj = parse_gbk(self.input_file)

        if obj.status == False:
            return False
        acs = obj.get_acs()

        for a in acs:
            features = obj.get_feature(a,self.feat_key)
            for seq in features:
                s = seq['dna_seq'].lower().replace('-','')
                for char in DNA_IUPAC_CHARS:
                    s = s.replace(char,"n")

                self.seq_data.append( {
                    'parent_id': a,
                    'locus_name':seq['gene_name'],
                    'seq_id': seq['gene_name'],
                    'dna_seq': s,
                    'dna_ambig_count': self.count_ambig_chars(seq['dna_seq'], DNA_AMBIG_CHARS),
                    'dna_hash': calc_md5([s])[0],
                    'dna_len': len(seq['dna_seq']),
                    'aa_seq': seq['aa_seq'],
                    'aa_hash': seq['aa_hash'],
                    'aa_len': len(seq['aa_seq']),

                } )
        return True

    def process_fasta(self):
        obj = parse_fasta(self.input_file)
        if obj.status == False:
            return
        ids = obj.get_seqids()
        for id in ids:
            features = obj.get_seq_by_id(id)
            seq = features['seq'].lower().replace('-','')
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

            self.seq_data.append({
                'parent_id': features['gene_name'],
                'locus_name': features['gene_name'],
                'seq_id': features['seq_id'],
                'dna_seq': dna_seq,
                'dna_ambig_count': self.count_ambig_chars(dna_seq, DNA_AMBIG_CHARS),
                'dna_hash': dna_hash,
                'dna_len': dna_len,
                'aa_seq': aa_seq,
                'aa_hash': aa_hash,
                'aa_len': aa_len,

            })

        return

    def process_seq_hash(self,sequences):
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
            self.seq_data.append({
                'parent_id': id,
                'locus_name': id,
                'seq_id': id,
                'dna_seq': dna_seq,
                'dna_hash': dna_hash,
                'dna_ambig_count':self.count_ambig_chars(dna_seq, DNA_AMBIG_CHARS),
                'dna_len': dna_len,
                'aa_seq': aa_seq,
                'aa_hash': aa_hash,
                'aa_len': aa_len,

            })

        return

    def count_ambig_chars(self,seq,chars):
        count = 0
        for char in chars:
            count+= seq.count(char)
        return count




class seq_store:
    stored_fields = ['parent_id','locus_name','seq_id','dna_hash','dna_len','aa_hash','aa_len','start_codon','stop_codon','count_internal_stop','dna_ambig_count']
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

    def __init__(self,sample_name,db_config_dict,metadata_dict,query_seq_records,blast_columns,filters={},stored_fields=[]):
        self.sample_name = sample_name
        self.record['query_data']['sample_name'] = sample_name
        self.add_db_config(db_config_dict)
        self.add_seq_data(query_seq_records)
        if len(stored_fields) > 0 :
            self.stored_fields = stored_fields
        self.add_db_metadata(metadata_dict)
        self.add_hit_cols(blast_columns)
        self.filters = filters
        self.record['query_data']['sample_name'] = self.sample_name


    def add_db_config(self,conf):
        self.record['db_info'] = conf

    def add_hit_cols(self,columns):
        self.record['query_hit_columns'] = columns

    def add_seq_data(self,query_seq_records):
        for idx in range(0, len(query_seq_records)):
            self.record['query_data']['query_seq_data'][idx] = {}
            for f in self.stored_fields:
                self.record['query_data']['query_seq_data'][idx][f] = ''
                if f in query_seq_records[idx]:
                    self.record['query_data']['query_seq_data'][idx][f] = query_seq_records[idx][f]

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






