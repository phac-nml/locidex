import hashlib
import json
import os
from collections import Counter

from Bio.Seq import Seq

from locidex.constants import NT_SUB, PROTEIN_ALPHA, DNA_ALPHA


def revcomp(s):
    """
    Helper method to create the reverse complement nucleotide sequence
    for the sequence passed.

    Parameters
    ----------
    s : str
        The nucleotide sequence to parse

    Returns
    -------
    str
        The reverse complement of the passed nucleotide sequence
    """
    return s.translate(NT_SUB)[::-1]

def calc_md5(list):
    """
    Method to encode the MD5 hash for the input string.

    Parameters
    ----------
    string : srt
    The string to compute the MD5 hash

    Returns
    -------
    hash
        The md5 hash generated
    """
    out = []
    for string in list:
        seq = str(string).encode()
        md5 = hashlib.md5()
        md5.update(seq)
        out.append(md5.hexdigest())
    return out

def translate_dna(dna_seq,trans_table=11):
    l = len(dna_seq)
    r = l % 3
    if r == 0:
        s = dna_seq
    else:
        s = ''.join(list(dna_seq)[:-r])
    return str(Seq(s).translate(table=trans_table))

def six_frame_translation(dna_seq,trans_table):
    fwd = []
    rev = []
    for i in range(0,3):
        fwd.append(translate_dna(dna_seq[i:], trans_table))
        rev.append(translate_dna(revcomp(dna_seq)[i:], trans_table))
    return (fwd, rev)

def guess_alphabet(seq):
    seq = seq.lower().replace('-','')
    l = len(seq)
    if l == 0:
        return None
    pr_inf_chars = set(PROTEIN_ALPHA) - set(DNA_ALPHA)
    counts = dict(Counter(seq))
    dtype = None
    if len(set(counts.keys()) & pr_inf_chars) > 0:
        dtype = 'protein'
    else:
        dna_bases = ['a','t','c','g','n']
        l = len(seq)
        c = 0
        for b in counts:
            if b not in dna_bases:
                continue
            c+= counts[b]
        if c / l >= 0.7:
            dtype = 'dna'
        else:
            dtype = 'protein'
    return dtype

def write_seq_list(seqs,output_file,format='json',seq_type='dna',seq_id_key='index'):
    with open(output_file, 'w') as oh:
        if format == 'json':
            json.dump(seqs, oh, indent=2)
        elif format == 'fasta':
            i = 0
            for record in seqs:
                if seq_id_key == 'index':
                    id = i
                else:
                    id = record[seq_id_key]
                s = record['dna_seq']
                if seq_type == 'aa':
                    s = record['aa_seq']
                oh.write(f">{id}\n{s}\n")

def write_seq_dict(data,output_file):
    with open(output_file, 'w') as oh:
        for id in data:
            oh.write(f">{id}\n{data[id]}\n")


def validate_dict_keys(data_dict,required_keys):
    for k in required_keys:
        if not k in data_dict:
            return [False, k]
    return [True, k ]

def get_file_length(self):
    return int(os.popen(f'wc -l {self.file_path}').read().split()[0])

def filter_hsps_df(df):
    return


