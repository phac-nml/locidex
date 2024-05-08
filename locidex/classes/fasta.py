from Bio import SeqIO
import gzip
from mimetypes import guess_type
from functools import partial
import os
from locidex.utils import calc_md5, slots
from dataclasses import dataclass

@dataclass
class Fasta:
    """
    """
    gene_name: str
    seq_id: str
    seq: str
    hash: str
    length: int
    __slots__ = slots(__annotations__)


class parse_fasta:

    def __init__(self, input_file,parse_def=False,seq_type=None,delim="|"):
        self.input_file = input_file
        if not os.path.isfile(self.input_file):
            raise FileNotFoundError("Input file: {} not found.".format(self.input_file))

        self.delim = delim
        self.seq_type = seq_type
        self.parse_def = parse_def
        self.seq_obj = self.parse_fasta()

    @staticmethod
    def normalize_sequence(fasta:str) -> str:
        """
        Remove INDELS and lower all characters in sequence
        """
        return fasta.lower().replace("-", "")

    def get_seqids(self):
        if self.seq_obj:
            return list(self.seq_obj.keys())
        raise AssertionError("No fasta file loaded.")

    def get_seq_by_id(self, fasta_id: str):
        if seq_data := self.seq_obj.get(fasta_id): #is not None and id in self.seq_obj:
            return seq_data
        raise KeyError("Missing sequence id: {}".format(fasta_id))

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
                data[id] = Fasta(gene_name=gene_name, seq_id=seq_id, seq=self.normalize_sequence(seq), hash=calc_md5([seq])[0], length=len(seq))
        return data
