import pandas as pd

from locidex.classes import run_command
import os

class blast_search:
    VALID_BLAST_METHODS = ['blastn','tblastn','blastp']
    BLAST_TABLE_COLS = []
    input_db_path = None
    input_query_path = None
    output_db_path = None
    output_results = None
    blast_params = {}
    blast_method = None
    create_db = False
    parse_seqids = False
    status = True
    messages = []

    def __init__(self,input_db_path,input_query_path,output_results,blast_params,blast_method,blast_columns,output_db_path=None,create_db=False,parse_seqids=False):
        self.input_query_path = input_query_path
        self.input_db_path = input_db_path
        self.output_db_path = output_db_path
        self.output_results = output_results
        self.blast_method = blast_method
        self.blast_params = blast_params
        self.create_db = create_db
        self.BLAST_TABLE_COLS = blast_columns
        self.parse_seqids = parse_seqids

        if self.output_db_path is None:
            self.output_db_path = input_db_path

        if not os.path.isfile(self.input_query_path):
            self.messages.append(f'Error {self.input_query_path} query fasta does not exist')
            self.status = False

        elif not create_db and not self.is_blast_db_valid():
            self.messages.append(f'Error {self.input_db_path} is not a valid blast db and creation of db is disabled')
            self.status = False
            return

        elif not blast_method in self.VALID_BLAST_METHODS:
            self.messages.append(f'Error {blast_method} is not a supported blast method: {self.blast_method}')
            self.status = False

        elif create_db and not self.is_blast_db_valid():
            (stdout,stderr) = self.makeblastdb()
            if not self.is_blast_db_valid():
                self.messages.append(f'Error {self.output_db_path} is not a valid blast db and creation of db failed')
                self.status = False


        self.messages.append(self.run_blast())



    def makeblastdb(self):
        dbtype = 'nucl'
        if self.blast_method == 'blastp':
            dbtype = 'prot'
        command = ['makeblastdb',
                   '-in', self.input_db_path,
                   '-dbtype', dbtype]
        if self.parse_seqids:
            command.append('-parse_seqids')
        if self.output_db_path != None:
            command +=['-out', self.output_db_path]

        return run_command(" ".join([str(x) for x in command]))

    def is_blast_db_valid(self):
        extensions = ['nsq', 'nin', 'nhr']
        for e in extensions:
            if not os.path.isfile(f'{self.input_db_path}.{e}'):
                extensions2 = ['pto', 'ptf', 'phr']
                for e2 in extensions2:
                    if not os.path.isfile(f'{self.input_db_path}.{e2}'):
                        return False

        return True

    def run_blast(self):
        command = [
            self.blast_method,
            '-query', self.input_query_path,
            '-db', self.input_db_path,
            '-out', self.output_results,
            '-outfmt', "'6 {}'".format(' '.join(self.BLAST_TABLE_COLS)),
        ]
        for p in self.blast_params:
            if p == 'parse_seqids':
                command.append(f'-{p}')
            else:
                command += [f'-{p}', f'{self.blast_params[p]}']
        return run_command(" ".join([str(x) for x in command]))


class parse_blast:
    BLAST_TABLE_COLS = []
    input_file = None
    df = None
    columns = []
    filter_options = {}
    status = True
    messages = []

    def __init__(self, input_file,blast_columns,filter_options):
        self.input_file = input_file
        self.filter_options = filter_options
        self.BLAST_TABLE_COLS = blast_columns

        if not os.path.isfile(self.input_file):
            self.messages.append(f'Error {self.input_file} does not exist')
            self.status = False
        self.read_hit_table()
        self.columns = self.df.columns.tolist()

        for id_col in ['qseqid','sseqid']:
            tp = {}
            if id_col in self.columns:
                tp[id_col] = 'object'
            self.df = self.df.astype(tp)
        for col_name in self.filter_options:
            if col_name in self.columns:
                min_value = self.filter_options[col_name]['min']
                max_value = self.filter_options[col_name]['max']
                include = self.filter_options[col_name]['include']
                self.filter_df(col_name, min_value, max_value,include)

    def read_hit_table(self):
        self.df = pd.read_csv(self.input_file,header=None,names=self.BLAST_TABLE_COLS,sep="\t",low_memory=False)


    def filter_df(self,col_name,min_value,max_value,include):
        if col_name not in self.columns:
            return False
        if min_value is not None:
            self.df = self.df[self.df[col_name] >= min_value]
        if max_value is not None:
            self.df = self.df[self.df[col_name] <= max_value]
        if include is not None:
            self.df = self.df[self.df[col_name].isin(include)]
        return True






