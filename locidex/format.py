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
from typing import List, Tuple
import logging
import errno
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pyrodigal import GeneFinder

from locidex.constants import FILE_TYPES, LocidexDBHeader, CharacterConstants
from locidex.utils import six_frame_translation, revcomp, calc_md5
from locidex.version import __version__

logger = logging.getLogger(__name__)
logging.basicConfig(filemode=sys.stderr, level=logging.INFO)

class locidex_format:

    input_type = None
    delim = '_'
    status = True
    __stop_codon = CharacterConstants.stop_codon

    # ? These two parameters below can probably be cleaned up
    __file_input = "file"
    __dir_input = "dir"

    @dataclass
    class FrameSelection:
        offset: int
        seq: str
        count_int_stops: int
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
        if self.input_type == self.__dir_input:
            self.process_dir()
        else:
            self.parse_fasta(self.input)

    def process_dir(self):
        files = self.get_dir_files(self.input)
        for f in files[self.__file_input]:
            for e in self.valid_ext:
                if e in f[1]:
                    self.gene_name = f[1].replace(f'{e}','')
                    self.parse_fasta(f[0])
                    break

    def set_input_type(self):
        if os.path.isfile(self.input):
            self.input_type = self.__file_input
        elif os.path.isdir(self.input):
            self.input_type = self.__dir_input
        else:
            logger.critical("Could not determine input type for: {}".format(self.input))
            raise AttributeError("Unknown input type could not be determined for: {}".format(self.input))

    def get_dir_files(self, input_dir):
        files = {self.__file_input: [], self.__dir_input: []}
        d = pathlib.Path(input_dir)
        for item in d.iterdir():
            type = self.__file_input
            if item.is_dir():
                type = self.__dir_input

            files[type].append([f"{item.resolve()}", os.path.basename(item)])
        return files

    def pick_frame(self, six_frame_translation) -> FrameSelection:
        """
        Reducing the complexity of this function now only checking if the allele is reverse complimented
        """

        reversed_frame_idx = 3 # all frames above this index are reverse complimented
        fwd_idx, rev_idx = 0, 3
        reverse_p = False

        rev = six_frame_translation[fwd_idx] # frame 1
        fwd = six_frame_translation[rev_idx] # frame 2

        fwd_stop_counts = fwd[:-1].count(self.__stop_codon)

        output_seq = fwd
        idx = fwd_idx
        min_int_stop = fwd_stop_counts

        if rev[-1] == self.__stop_codon and (stop_counts := rev[:-1].count(self.__stop_codon)) == 1 and fwd_stop_counts > 0:
            idx = rev_idx
            output_seq = rev
            min_int_stop = stop_counts

        offset = idx % reversed_frame_idx # gives offset for both revcomp and seq
        if idx >= reversed_frame_idx:
            reverse_p = True

        return self.FrameSelection(offset=offset, seq=output_seq.lower(), count_int_stops=min_int_stop, revcomp=reverse_p)


    def parse_fasta(self, input_file):
        encoding = guess_type(input_file)
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(input_file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                id = str(record.id)
                if self.input_type == self.__file_input:
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
        logger.critical(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        raise FileExistsError(errno.EEXIST, os.strerror(errno.EEXIST), str(outdir))


    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    if not os.path.isdir(input) and not os.path.isfile(input):
        logger.critical(f'Error {input} does not exist as a file or directory')
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), str(input))
    logger.info("Beginning format operation.")
    obj = locidex_format(input=input,header=LocidexDBHeader._fields,is_protein=is_coding,min_len_frac=min_len_frac,max_len_frac=max_len_frac, min_ident_perc=min_ident,
            min_cov_perc=min_match_cov,trans_table=trans_table,valid_ext=FILE_TYPES['fasta'])
    logger.info("Finished format.")
    run_data['result_file'] = os.path.join(outdir,"locidex.txt")
    pd.DataFrame.from_dict(obj.data,orient='index').to_csv(run_data['result_file'],sep="\t",index=False,header=True)

    run_data['analysis_end_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(os.path.join(outdir,"results.json"),"w") as oh:
        oh.write(json.dumps(run_data,indent=4))


# call main function
if __name__ == '__main__':
    run()



