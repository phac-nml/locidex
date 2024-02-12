import json
from datetime import datetime
import pandas as pd
import os, sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from locidex.version import __version__
from locidex.constants import FORMAT_RUN_DATA
from locidex.classes import run_command

class locidex_build:
    input_file = None
    outdir = None
    blast_dir = None
    df = None
    meta = {}
    status = True
    messages = []


    def __init__(self, input_file, outdir,seq_columns=['dna_seq','aa_seq'],force=False):
        self.input_file = input_file
        self.outdir = outdir
        self.force = force
        self.status = self.init_dir(self.outdir)
        if not self.status:
            return

        self.status = self.blast_dir = os.path.join(outdir,'blast')
        self.init_dir(self.blast_dir)
        if not self.status:
            return

        self.df = self.read_data( self.input_file)

    def init_dir(self,d):
        if os.path.isdir(d) and not self.force:
            self.messages.append(f'Error {d} exists, if you would like to overwrite, then specify --force')
            return False

        if not os.path.isdir(self.d):
            os.makedirs(self.d, 0o755)

        if not os.path.isdir(self.outdir):
            os.makedirs(self.d, 0o755)

        return True

    def read_data(self,f):
        if os.path.isfile(f) and os.path.getsize(f) > 0:
            return pd.read_csv(f,sep="\t",header=0)
        else:
            return pd.DataFrame()

    def set_coding(self,col_name):
        s = True
        if not col_name in self.df:
            return False
        values = self.df[col_name].astype('str').unique

    def write_seq_file(self,f,col_name):
        oh = open(f, "w")
        data = dict(zip(self.df.index.tolist(), self.df[col_name]))
        for id in data:
            oh.write(f">{id}\n{data[id]}\n")
        oh.close()

    def makeblastdb(self):
        dbtype = 'nucl'
        if self.blast_method == 'blastp':
            dbtype = 'prot'
        command = ['makeblastdb',
                   '-in', self.input_query_path,
                   '-dbtype', dbtype]
        if self.parse_seqids:
            command.append('-parse_seqids')
        if self.output_db_path != None:
            command.append('-out',self.output_db_path)

        return run_command(command)













def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Locidex: Format an existing allele database into Locidex build tsv format",
        formatter_class=CustomFormatter)
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

    return parser.parse_args()


def run():
    cmd_args = parse_args()
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

    run_data = FORMAT_RUN_DATA
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)

    run_data['result_file'] = os.path.join(outdir)

    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    with open(os.path.join(outdir,"results.json"),"w") as out:
        json.dump(run_data,indent=4)






# call main function
if __name__ == '__main__':
    run()



