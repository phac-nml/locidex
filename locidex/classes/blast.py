import pandas as pd

from locidex.classes import run_command
from locidex.utils import slots
from dataclasses import dataclass
from typing import Optional
import os


@dataclass
class FilterOptions:
    min: Optional[float]
    max: Optional[float] 
    include: Optional[bool] 
    __slots__ = slots(__annotations__)

class blast_search:
    VALID_BLAST_METHODS = ['blastn','tblastn','blastp']

    def __init__(self,input_db_path,input_query_path,output_results,blast_params,blast_method,blast_columns,output_db_path=None,parse_seqids=False):
        self.input_query_path = input_query_path
        self.input_db_path = input_db_path
        self.output_db_path = output_db_path
        self.output_results = output_results
        self.blast_method = blast_method
        self.blast_params = blast_params
        #self.create_db = create_db
        self.BLAST_TABLE_COLS = blast_columns
        self.parse_seqids = parse_seqids

        if self.output_db_path is None:
            self.output_db_path = input_db_path

        if not os.path.isfile(self.input_query_path):
            raise ValueError(f'Error {self.input_query_path} query fasta does not exist')
        elif not self.is_blast_db_valid():
            raise ValueError(f'Error {self.input_db_path} is not a valid blast db')
        elif not blast_method in self.VALID_BLAST_METHODS:
            raise ValueError(f'Error {blast_method} is not a supported blast method: {self.blast_method}')

        #elif create_db and not self.is_blast_db_valid():
        #    (stdout,stderr) = self.makeblastdb()
        #    print(stdout, stderr)
        #    if not self.is_blast_db_valid():
        #        raise ValueError(f'Error {self.output_db_path} is not a valid blast db and creation of db failed')
        # TODO add logging for blast messages
        (stdout, stderr) = self.run_blast()
        print(stdout, stderr)


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
            print(self.output_db_path, e)
            if not os.path.isfile(f'{self.output_db_path}.{e}'):
                extensions2 = ['pto', 'ptf', 'phr']
                for e2 in extensions2:
                    if not os.path.isfile(f'{self.output_db_path}.{e2}'):
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
    input_file = None
    columns = []
    filter_options = {}

    def __init__(self, input_file,blast_columns,filter_options):
        self.input_file = input_file
        self.filter_options = filter_options
        self.BLAST_TABLE_COLS = blast_columns

        if not os.path.isfile(self.input_file):
            raise FileNotFoundError(f'Error {self.input_file} does not exist')
            
        self.df = self.read_hit_table()
        print(self.df)
        self.columns = self.df.columns.tolist()

        for id_col in ['qseqid','sseqid']:
            tp = {}
            if id_col in self.columns:
                tp[id_col] = 'object'
            self.df = self.df.astype(tp)
        for col_name in self.filter_options:
            if col_name in self.columns:
                #min_value = self.filter_options[col_name]['min']
                #max_value = self.filter_options[col_name]['max']
                #include = self.filter_options[col_name]['include']
                min_value = self.filter_options[col_name].min
                max_value = self.filter_options[col_name].max
                include = self.filter_options[col_name].include
                self.filter_df(col_name, min_value, max_value,include)

    def read_hit_table(self):
        return pd.read_csv(self.input_file,header=None,names=self.BLAST_TABLE_COLS,sep="\t",low_memory=False)


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






