import json
from datetime import datetime
import pandas as pd
import os, sys
from pathlib import Path
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from locidex.version import __version__
from locidex.constants import DBFiles
from locidex.classes import run_command
from locidex.constants import DBConfig, MetadataFields, raise_file_not_found_e
from locidex.classes.blast import BlastMakeDB
from locidex.manifest import DBData
import getpass
import errno
import logging
import sys

logger = logging.getLogger(__name__)
logging.basicConfig(filemode=sys.stderr, level=logging.INFO)

class locidex_build:


    def __init__(self, input_file: os.PathLike, outdir: os.PathLike, config: DBConfig, 
                seq_columns={DBData.nucleotide_name():'dna_seq',DBData.protein_name():'aa_seq'},force=False,parse_seqids=False,translation_table: int = 11):
        self.input_file = input_file
        self.translation_table = translation_table
        self.outdir = outdir
        self.force = force
        self.status = self.init_dir(self.outdir)
        self.parse_seqids = parse_seqids

        self.is_dna = False
        self.is_protein = False
        self.blast_dir = outdir.joinpath("blast")

        self.init_dir(self.blast_dir)

        self.df = self.read_data(self.input_file)
        self.config = config
        self.config.db_num_seqs = len(self.df)

        for t in seq_columns:
            col_name = seq_columns[t]
            s = self.is_seqtype_present(col_name)
            if s:
                outfile  = self.blast_dir.joinpath(t)
                if t == DBData.nucleotide_name():
                    self.is_dna = True
                    self.config.nucleotide_db_name = t
                    blast_method = DBData.nucleotide_db_type()
                elif t == DBData.protein_name():
                    self.is_protein  = True
                    self.config.protein_db_name = t
                    blast_method = DBData.protein_db_type()
                self.init_dir(outfile)
                out_fasta = outfile.joinpath("{}.fasta".format(t))
                self.write_seq_file(out_fasta, col_name)
                creating_db = BlastMakeDB(input_file=out_fasta, 
                                        db_type=blast_method,
                                        parse_seqids=self.parse_seqids,
                                        output_db_path=out_fasta)
                creating_db.makeblastdb()
        self.config.is_nucl = self.is_dna
        self.config.is_prot = self.is_protein
        self.get_metadata(self.df,columns_to_exclude=list(seq_columns.values()))

    def init_dir(self, d: Path):
        logger.info("Checking if directory {} exists.".format(d))
        try:
            d.mkdir(parents=True, exist_ok=self.force, mode=0o755)
        except FileExistsError:
            logger.critical("Database file {} already exists. To overwrite please run with --force".format(d))
            raise FileExistsError(errno.EEXIST, os.strerror(errno.EEXIST), str(d))
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

    def get_metadata(self,df,columns_to_exclude=['dna_seq','aa_seq']):
        columns_to_exclude = set(columns_to_exclude)
        columns_to_include = []
        for col in list(df.columns):
            if col not in columns_to_exclude:
                columns_to_include.append(col)
        subset = df[columns_to_include]
        self.meta = {
            "info": MetadataFields(
                num_seqs=len(df),
                is_cds=True,
                trans_table=self.translation_table,
                dna_min_len=min(df['dna_min_len'].tolist()),
                dna_max_len=max(df['dna_min_len'].tolist()),
                dna_min_ident=min(df['dna_min_ident'].tolist()),
                aa_min_len=min(df['aa_min_len'].tolist()),
                aa_max_len= max(df['aa_min_len'].tolist()),
                aa_min_ident= min(df['aa_min_ident'].tolist())).to_dict(),
            "meta": subset.to_dict(orient='index')
        }

def add_args(parser=None):
    translation_table_def = 11
    default_version = "1.0.0"
    if parser is None:
        parser = ArgumentParser(
            description="Locidex: Build a locidex database",)
    parser.add_argument('-i','--input_file', type=str, required=True,help='Input tsv formated for locidex')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output directory to put results')
    parser.add_argument('-n', '--name', type=str, required=True, help='Database name',default='Locidex Database')
    parser.add_argument('-a', '--author', type=str, required=False, help='Author Name for Locidex Database. Default: {}'.format(getpass.getuser()),default=getpass.getuser())
    parser.add_argument('-c', '--db_ver', type=str, required=False, help='Version code for locidex db: {}'.format(default_version),
                        default=default_version)
    parser.add_argument('-e', '--db_desc',type=str, required=False, help='Version code for locidex db',
                        default='')
    parser.add_argument("-t", "--translation_table", type=int, required=False, help="Translation table to use: {}".format(translation_table_def), default=translation_table_def)
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    return parser
    



def run(cmd_args=None):
    if cmd_args is None:
        parser = add_args()
        cmd_args = parser.parse_args()
    
    input_file = cmd_args.input_file
    outdir = cmd_args.outdir
    force = cmd_args.force


    config = DBConfig(
        db_name=cmd_args.name,
        db_version =cmd_args.db_ver,
        db_desc=cmd_args.db_desc,
        db_author=cmd_args.author,
        db_date=datetime.now().strftime("%Y/%d/%m"),
    )

    run_params = vars(cmd_args)
    run_params = run_params | config.to_dict()
    run_data = dict()
    run_data['analysis_start_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)

    if not os.path.isfile(input_file):
        logger.critical(f'Error {input_file} does not exist, please check path and try again')
        raise_file_not_found_e(input_file, logger)

    obj = locidex_build(Path(input_file), Path(outdir),config=config,seq_columns={'nucleotide':'dna_seq','protein':'aa_seq'},force=force)

    if obj.status == False:
        logger.critical(f'Error something went wrong building the db.')
        raise ValueError("Something went wrong building db.")

    run_data['analysis_end_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(os.path.join(outdir,DBFiles.config_file),"w") as oh:
        oh.write(json.dumps(obj.config.to_dict(),indent=4))

    with open(os.path.join(outdir, DBFiles.meta_file),"w") as oh:
        oh.write(json.dumps(obj.meta,indent=4))


    with open(os.path.join(outdir, DBFiles.results_file),"w") as oh:
        oh.write(json.dumps(run_data,indent=4))


# call main function
if __name__ == '__main__':
    run()



