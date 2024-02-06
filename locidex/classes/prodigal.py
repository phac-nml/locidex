from Bio import SeqIO
import os
import gzip
from mimetypes import guess_type
from functools import partial
import pyrodigal

class gene_prediction:
    file = None
    sequences = {}
    status = True
    messages = []

    def __init__(self,file):
        self.file = file
        if not os.path.isfile(self.file):
            self.messages.append(f'Error {self.file} does not exist')
            self.status = False


    def predict(self):
        encoding = guess_type(self.file)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(self.file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                orf_finder = pyrodigal.GeneFinder(meta=True)
                for i, pred in enumerate(orf_finder.find_genes(bytes(record.seq))):
                    id = f">{record.id}_{i + 1}"
                    self.sequences[id] = pred.sequence()

