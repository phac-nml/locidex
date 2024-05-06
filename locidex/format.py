import gzip
import json
import os
import pathlib
import sys
from argparse import ArgumentParser
from datetime import datetime
from functools import partial
from mimetypes import guess_type
from dataclasses import dataclass

import pandas as pd
from Bio import SeqIO

from locidex.constants import FILE_TYPES, LocidexDBHeader
from locidex.utils import six_frame_translation, revcomp, calc_md5
from locidex.version import __version__

class locidex_format:

    input_type = None
    delim = '_'
    status = True

    @dataclass
    class FrameSelection:
        offset: int
        seq: str
        count_int_stops: int
        frame: int
        revcomp: bool

    def __init__(self,input,header,is_protein=False,delim="_",trans_table=11,
                min_len_frac=0.7,max_len_frac=1.3,min_cov_perc=80.0,min_ident_perc=80.0,valid_ext=None):
        self.input = input
        self.seq_idx = 0
        self.gene_name = None
        self.header = header
        self.delim = delim
        self.translation_table = trans_table
        self.min_len_frac = min_len_frac
        self.max_len_frac = max_len_frac
        self.min_cov_perc = min_cov_perc
        self.min_ident_perc = min_ident_perc
        self.is_protein_coding = is_protein
        self.data = dict()
        self.valid_ext = valid_ext

        if valid_ext is not None:
            if isinstance(valid_ext,list):
                self.valid_ext = valid_ext
            else:
                self.valid_ext = [valid_ext]

        self.set_input_type()
        if self.input_type == 'dir':
            self.process_dir()
        else:
            self.parse_fasta(self.input)

    def process_dir(self):
        files = self.get_dir_files(self.input)
        for f in files['file']:
            for e in self.valid_ext:
                if e in f[1]:
                    self.gene_name = f[1].replace(f'.{e}','')
                    self.parse_fasta(f[0])
                    break

    def set_input_type(self):
        if os.path.isfile(self.input):
            self.input_type = 'file'
        elif os.path.isdir(self.input):
            self.input_type = 'dir'
        else:
            raise AttributeError("Unknown input type could not be determined for: {}".format(self.input))

    def get_dir_files(self, input_dir):
        files = {'file': [], 'dir': []}
        d = pathlib.Path(input_dir)
        for item in d.iterdir():
            type = 'file'
            if item.is_dir():
                type = 'dir'

            files[type].append([f"{item.resolve()}", os.path.basename(item)])
        return files

    def pick_frame(self, six_frame_translation) -> FrameSelection:
        count_internal_stops = []
        terminal_stop_codon_present = []
        for i in range(0, len(six_frame_translation)):
            for k in range(0, len(six_frame_translation[i])):
                count_internal_stops.append(six_frame_translation[i][k][:-1].count('*'))
                terminal_stop_codon_present.append(six_frame_translation[i][k][-1] == '*')

        min_int_stop = min(count_internal_stops)
        idx = count_internal_stops.index(min_int_stop)
        i = 0
        k = 0
        offset = 0
        r = False
        s = 1

        if min_int_stop == 0 and (terminal_stop_codon_present[idx] == 1):
            s = idx + 1
            if s == 1 or s == 4:
                offset = 0
            elif s == 2 or s == 5:
                offset = 1
            else:
                offset = 2
            if idx > 2:
                r = True
            k = offset
        elif min_int_stop == 0 and max(terminal_stop_codon_present) == 1:
            best_idx = [0, min_int_stop, False]
            for idx, value in enumerate(count_internal_stops):
                if value != min_int_stop:
                    continue
                if best_idx[2] == False and terminal_stop_codon_present[idx]:
                    best_idx = [idx, min_int_stop, terminal_stop_codon_present[idx]]
            s = best_idx[0] + 1
            if s == 1 or s == 4:
                offset = 0
            elif s == 2 or s == 5:
                offset = 1
            else:
                offset = 2
            if idx > 2:
                r = True
            k = offset

        if idx > 2:
            i = 1

        seq = six_frame_translation[i][k].lower()
        return self.FrameSelection(offset=offset, frame=s, seq=seq, count_int_stops=min_int_stop, revcomp=r)


    def parse_fasta(self, input_file):
        encoding = guess_type(input_file)
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(input_file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                id = str(record.id)
                if self.input_type == 'file':
                    gene_name = "_".join(id.split(self.delim)[:-1])
                else:
                    gene_name = self.gene_name
                dna_seq = str(record.seq).lower().replace('-','')
                if self.is_protein_coding:
                    t = self.pick_frame(six_frame_translation(dna_seq, trans_table=self.translation_table))
                    aa_seq = t.seq
                    dna_seq = dna_seq[t.offset:]
                    if t.revcomp:
                        dna_seq = revcomp(dna_seq)

                dna_len = len(dna_seq)
                aa_len = len(aa_seq)

                aa_encoding_p = lambda x: x if self.is_protein_coding else None
                row = LocidexDBHeader(
                    seq_id=self.seq_idx,
                    locus_name=gene_name,
                    locus_name_alt=id,
                    locus_product='',
                    locus_description='',
                    locus_uid=id.split(self.delim)[-1],
                    dna_seq=dna_seq,
                    dna_seq_len=dna_len,
                    dna_seq_hash=calc_md5([dna_seq])[0],
                    dna_ambig_count=dna_seq.count('n'),
                    aa_seq = aa_encoding_p(aa_seq),
                    aa_seq_len=aa_encoding_p(len(aa_seq)),
                    aa_seq_hash= aa_encoding_p(calc_md5([aa_seq])[0]),
                    aa_min_ident=aa_encoding_p(self.min_ident_perc * (self.min_ident_perc /100)),
                    aa_min_len=aa_encoding_p(aa_len * self.min_len_frac),
                    aa_max_len=aa_encoding_p( aa_len * self.max_len_frac),
                    min_aa_match_cov=aa_encoding_p(self.min_cov_perc),
                    count_int_stops=aa_encoding_p(t.count_int_stops),
                    dna_min_len=dna_len*self.min_len_frac,
                    dna_max_len=dna_len*self.max_len_frac,
                    dna_min_ident=self.min_ident_perc,
                    min_dna_match_cov=self.min_cov_perc
                )

                self.data[self.seq_idx] = row
                self.seq_idx += 1


def add_args(parser=None):
    if parser is None:
        parser = ArgumentParser(
            description="Locidex: Format sequences for a database.",)
    parser.add_argument('-i','--input', type=str, required=True,help='Input directory of fasta files or input fasta')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output directory to put results')
    parser.add_argument('--min_len_frac', type=int, required=False, help='Used to calculate individual sequence minimum acceptable length (0 - 1)',
                        default=0.7)
    parser.add_argument('--max_len_frac', type=int, required=False, help='Used to calculate individual sequence maximimum acceptable length (1 - n)',
                        default=1.3)
    parser.add_argument('--min_ident', type=float, required=False, help='Global minumum percent identity required for match',
                        default=80.0)
    parser.add_argument('--min_match_cov', type=float, required=False, help='Global minumum percent hit coverage identity required for match',
                        default=80.0)
    parser.add_argument('--translation_table', type=int, required=False,
                        help='output directory', default=11)
    parser.add_argument('-n', '--not_coding', required=False, help='Skip translation',
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    return parser

def run(cmd_args=None):
    if cmd_args is None:
        parser = add_args()
        cmd_args = parser.parse_args()

    input = cmd_args.input
    outdir = cmd_args.outdir
    min_len_frac = cmd_args.min_len_frac
    max_len_frac = cmd_args.max_len_frac
    min_ident = cmd_args.min_ident
    min_match_cov = cmd_args.min_match_cov
    trans_table = cmd_args.translation_table
    force = cmd_args.force

    is_coding = True
    if cmd_args.not_coding:
        is_coding = False

    run_data = dict()
    run_data['analysis_start_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)

    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    if not os.path.isdir(input) and not os.path.isfile(input):
        print(f'Error {input} does not exist as a file or directory')
        sys.exit()

    obj = locidex_format(input=input,header=LocidexDBHeader._fields,is_protein=is_coding,min_len_frac=min_len_frac,max_len_frac=max_len_frac, min_ident_perc=min_ident,
            min_cov_perc=min_match_cov,trans_table=trans_table,valid_ext=FILE_TYPES['fasta'])

    run_data['result_file'] = os.path.join(outdir,"locidex.txt")
    pd.DataFrame.from_dict(obj.data,orient='index').to_csv(run_data['result_file'],sep="\t",index=False,header=True)

    run_data['analysis_end_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(os.path.join(outdir,"results.json"),"w") as oh:
        oh.write(json.dumps(run_data,indent=4))


# call main function
if __name__ == '__main__':
    run()



