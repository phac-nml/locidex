import os
from locidex.classes.gbk import parse_gbk
from locidex.classes.fasta import parse_fasta
from locidex.utils import guess_alphabet, calc_md5, six_frame_translation

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

    def __init__(self,input_file,file_type,feat_key='CDS',translation_table=11):
        self.input_file = input_file
        self.file_type = file_type
        self.translation_table = translation_table
        self.feat_key = feat_key

        if not os.path.isfile(self.input_file):
            self.messages.append(f'Error {self.input_file} does not exist')
            self.status = False
            return self.status

        if file_type == 'genbank':
            self.status = self.process_gbk()
        elif file_type == 'fasta':
            self.status = self.process_fasta()
        elif file_type == 'gff':
            self.status = False
        elif file_type == 'gtf':
            self.status = False

        if self.status:
            self.add_codon_data()

        return self.status


    def add_codon_data(self):
        for record in self.seq_data:
            if record['aa_len'] == 0:
                continue
            dna_len = record['dna_seq']
            dna_seq = record['dna_seq']
            start_codon = ''
            stop_codon = ''
            count_internal_stop = 0
            if dna_len >= 6:
                start_codon = dna_seq[0:3]
                stop_codon = dna_seq[-3:]
                count_internal_stop = record['aa_seq'].count('*')
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
                self.seq_data.append( {
                    'parent_id': a,
                    'locus_name':seq['gene_name'],
                    'seq_id': seq['gene_name'],
                    'dna_seq': seq['dna_seq'],
                    'dna_hash': seq['dna_hash'],
                    'dna_len': len(seq['dna_seq']),
                    'aa_seq': seq['aa_seq'],
                    'aa_hash': seq['aa_hash'],
                    'aa_len': len(seq['aa_seq']),

                } )
        return True

    def process_fasta(self):
        obj = parse_fasta(self.input_file)
        if obj.status == False:
            return False
        ids = obj.get_seqids()
        for id in ids:
            features = obj.get_feature(id, self.feat_key)
            for seq in features:
                seq = seq['dna_seq']
                dtype = guess_alphabet(seq)
                dna_seq = ''
                dna_hash = ''
                dna_len = 0
                aa_seq = ''
                aa_hash = ''
                aa_len = 0
                if dtype == 'dna':
                    dna_seq = seq
                    dna_hash = calc_md5([seq])[0]
                    dna_len = len(seq)
                    aa_seq = six_frame_translation(dna_seq,self.translation_table)
                    aa_hash = calc_md5([aa_seq])[0]
                    aa_len = len(aa_seq)
                else:
                    aa_seq = six_frame_translation(dna_seq,self.translation_table)[0][0]
                    aa_hash = calc_md5([aa_seq])[0]
                    aa_len = len(aa_seq)


                self.seq_data.append({
                    'parent_id': seq['gene_name'],
                    'locus_name': seq['gene_name'],
                    'seq_id': seq['seq_id'],
                    'dna_seq': dna_seq,
                    'dna_hash': dna_hash,
                    'dna_len': dna_len,
                    'aa_seq': aa_seq,
                    'aa_hash': aa_hash,
                    'aa_len': aa_len,

                })

        return True





