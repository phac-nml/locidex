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
    num_seqs = 0
    status = True
    messages = []


    def __init__(self, input_file, outdir,config={},seq_columns={'nucleotide':'dna_seq','protein':'aa_seq'},force=False,parse_seqids=False):
        self.input_file = input_file
        self.outdir = outdir
        self.force = force
        self.status = self.init_dir(self.outdir)
        self.parse_seqids = parse_seqids
        if not self.status:
            return

        self.status = self.blast_dir = os.path.join(outdir,'blast')
        self.init_dir(self.blast_dir)
        if not self.status:
            return

        self.df = self.read_data( self.input_file)
        self.config = config
        self.config["db_num_seqs"] = len(self.df)

        for t in seq_columns:
            col_name = seq_columns[t]
            s = self.is_seqtype_present(col_name)
            if s:
                outfile  = os.path.join(self.blast_dir, t)
                if t == 'nucleotide':
                    self.is_dna = True
                    self.config["nucleotide_db_name"] = t
                    blast_method = 'nucl'
                elif t == 'protein':

                    self.is_protein  = True
                    self.config["protein_db_name"] = t
                    blast_method = 'prot'
                self.create_seq_db(t, col_name, outfile, blast_method)
        self.config["is_nucl"] = self.is_dna
        self.config["is_prot"] = self.is_protein
        self.get_metadata(self.df,columns_to_exclude=list(seq_columns.values()))

    def create_seq_db(self,stype,col_name,outfile,blast_method='nucl'):
        self.init_dir(outfile)
        f = os.path.join(outfile,"{}.fasta".format(stype))
        self.write_seq_file(f,col_name)
        (stdout,stderr) = self.makeblastdb(fasta_file=f,outfile=os.path.join(outfile,stype),blast_method=blast_method)
        self.messages.append(stderr)

    def init_dir(self,d):
        if os.path.isdir(d) and not self.force:
            self.messages.append(f'Error {d} exists, if you would like to overwrite, then specify --force')
            return False

        if not os.path.isdir(d):
            os.makedirs(d, 0o755)

        if not os.path.isdir(d):
            os.makedirs(d, 0o755)

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
        values = list(self.df[col_name].dropna().astype('str').unique())
        if len(values) == 0:
            return False
        return True


    def write_seq_file(self,f,col_name):
        oh = open(f, "w")
        data = dict(zip(self.df.index.tolist(), self.df[col_name]))
        for id in data:
            oh.write(f">{id}\n{data[id]}\n")
        oh.close()

    def makeblastdb(self,fasta_file,outfile,blast_method):
        dbtype = 'nucl'
        if blast_method == 'prot':
            dbtype = 'prot'
        command = ['makeblastdb',
                   '-in', fasta_file,
                   '-dbtype', dbtype]
        if self.parse_seqids:
            command.append('-parse_seqids')
        if outfile != None:
            command+= ['-out',outfile]
        command = " ".join([str(x) for x in command])

        return run_command(command)

    def get_metadata(self,df,columns_to_exclude=['dna_seq','aa_seq']):
        columns_to_exclude = set(columns_to_exclude)
        columns_to_include = []
        for col in list(df.columns):
            if col not in columns_to_exclude:
                columns_to_include.append(col)
        subset = df[columns_to_include]
        self.meta = {
            "info": {
                "num_seqs": len(df),
                "is_cds": "True",
                "trans_table": 11,
                "dna_min_len": min(df['dna_min_len'].tolist()),
                "dna_max_len": max(df['dna_min_len'].tolist()),
                "dna_min_ident": min(df['dna_min_ident'].tolist()),
                "aa_min_len": min(df['aa_min_len'].tolist()),
                "aa_max_len": max(df['aa_min_len'].tolist()),
                "aa_min_ident": min(df['aa_min_ident'].tolist()),
            },
            "meta": subset.to_dict(orient='index')
        }













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
    parser.add_argument('-e', '--db_desc',type=str, required=False, help='Version code for locidex db',
                        default='')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')

    return parser.parse_args()


def run():
    cmd_args = parse_args()
    input_file = cmd_args.input_file
    outdir = cmd_args.outdir
    force = cmd_args.force
    run_data = FORMAT_RUN_DATA
    run_data['analysis_start_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
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
    if obj.status == False:
        print(f'Error something went wrong building the db, check error messages {obj.messages}')
        sys.exit()

    run_data['analysis_end_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(os.path.join(outdir,"config.json"),"w") as oh:
        oh.write(json.dumps(obj.config,indent=4))

    with open(os.path.join(outdir,"meta.json"),"w") as oh:
        oh.write(json.dumps(obj.meta,indent=4))


    with open(os.path.join(outdir,"results.json"),"w") as oh:
        oh.write(json.dumps(run_data,indent=4))






# call main function
if __name__ == '__main__':
    run()



