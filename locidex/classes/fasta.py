from Bio import SeqIO
import gzip
from mimetypes import guess_type
from functools import partial
import os
from locidex.utils import calc_md5


class parse_fasta:
    input_file = None
    seq_obj = None
    status = True
    messages = []

    def __init__(self, input_file,parse_def=False,seq_type=None,delim="|"):
        self.input_file = input_file
        if not os.path.isfile(self.input_file):
            self.messages.append(f'Error {self.input_file} does not exist')
            self.status = False

        self.delim = delim
        self.seq_type = seq_type
        self.parse_def = parse_def
        self.seq_obj = self.parse_fasta()


        return

    def get_seqids(self):
        if self.seq_obj is not None:
            return list(self.seq_obj.keys())
        else:
            return []

    def get_seq_by_id(self, id):
        if self.seq_obj is not None and id in self.seq_obj:
            return self.seq_obj[id]
        else:
            return {}


    def parse_fasta(self):
        encoding = guess_type(self.input_file)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        data = {}
        with _open(self.input_file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                id = str(record.id)
                seq = str(record.seq)
                gene_name = id
                seq_id = id
                if self.parse_def:
                    h = id.split(self.delim)
                    if len(h) > 1:
                        gene_name = h[0]
                        seq_id = h[1]
                data[id] = {'gene_name': gene_name, 'seq_id':seq_id,'seq': seq, 'hash':calc_md5([seq])[0],'length':len(seq)}
        return data
