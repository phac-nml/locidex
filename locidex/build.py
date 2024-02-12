import json
from datetime import datetime
import pandas as pd
import os, sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from locidex.version import __version__
from locidex.constants import FORMAT_RUN_DATA, DB_CONFIG_FIELDS
from locidex.classes import run_command

class locidex_build:
    input_file = None
    outdir = None
    blast_dir = None
    is_dna = False
    is_protein = False
    df = None
    config = {}
    meta = {}
    status = True
    messages = []


    def __init__(self, input_file, outdir,config={},seq_columns={'nucleotide':'dna_seq','protein':'aa_seq'},force=False):
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
        self.config["db_num_seqs"] = len(self.df)
        for t in seq_columns:
            col_name = seq_columns[t]
            s = self.is_seqtype_present(col_name)
            if s:
                self.create_seq_db(t, col_name)
                if t == 'nucleotide':
                    self.is_dna = True
                    self.config["nucleotide_db_name"] = t
                elif t == 'protien':
                    self.is_protein  = True
                    self.config["protein_db_name"] = t
        self.config["is_nucl"] = self.is_dna
        self.config["is_prot"] = self.is_protein

    def create_seq_db(self,stype,col_name):
        d = os.path.join(self.blast_dir,stype)
        self.init_dir(d)
        f = os.path.join(d,"{}.fasta".format(stype))
        self.write_seq_file(f,col_name)


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

    def is_seqtype_present(self,col_name):
        s = True
        if not col_name in self.df:
            return False
        values = self.df[col_name].dropna().astype('str').unique().to_list()
        if len(values) == 0:
            return False
        return True


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
        description="Locidex: Build a locidex database",
        formatter_class=CustomFormatter)
    parser.add_argument('-i','--input_file', type=str, required=True,help='Input tsv formated for locidex')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output directory to put results')
    parser.add_argument('-n', '--name', type=str, required=False, help='DB name',default='Locidex Database')
    parser.add_argument('-a', '--author', type=str, required=False, help='Author Name for Locidex Database',default='')
    parser.add_argument('-d', '--date', type=str, required=False, help='Creation date for Locidex Database',
                        default='')
    parser.add_argument('-c', '--db_ver', type=str, required=False, help='Version code for locidex db',
                        default='1.0.0')
    parser.add_argument('-c', '--db_desc',type=str, required=False, help='Version code for locidex db',
                        default='')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')

    return parser.parse_args()


def run():
    cmd_args = parse_args()
    input_file = cmd_args.input
    outdir = cmd_args.outdir
    force = cmd_args.force
    run_data = FORMAT_RUN_DATA
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)

    config = {}
    for f in DB_CONFIG_FIELDS:
        config[f] = ''

    config["db_name"] = cmd_args.name
    config["db_version"] = cmd_args.db_ver
    config["db_desc"] = cmd_args.db_desc
    config["db_author"] = cmd_args.author
    if cmd_args.date == '':
        config["db_date"] = datetime.now().strftime("%d/%m/%Y")

    if not os.path.isfile(input_file):
        print(f'Error {input_file} does not exist, please check path and try again')
        sys.exit()

    run_data['result_file'] = os.path.join(outdir)
    obj = locidex_build(input_file, outdir,config=config,seq_columns={'nucleotide':'dna_seq','protein':'aa_seq'},force=force)
    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    with open(os.path.join(outdir,"results.json"),"w") as out:
        json.dump(run_data,indent=4)






# call main function
if __name__ == '__main__':
    run()



