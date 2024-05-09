"""
Blast module refactored
"""

import pandas as pd

from locidex.classes import run_command
from locidex.utils import slots
from locidex.manifest import DBData
from locidex.constants import BlastCommands
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict
import logging
import sys

logger = logging.getLogger(__name__)
logging.basicConfig(filemode=sys.stderr, level=logging.INFO)

@dataclass
class FilterOptions:
    min: Optional[float]
    max: Optional[float] 
    include: Optional[bool] 
    __slots__ = slots(__annotations__)


class BlastMakeDB:
    """
    Create a blast database
    """

    def __init__(self, input_file: Path, db_type: str, parse_seqids: bool, output_db_path: Optional[Path]):
        self.input_file = input_file
        self.db_type = db_type
        self.parse_seqids = parse_seqids
        self.output_db_path = output_db_path

    def makeblastdb(self):
        command = ['makeblastdb', '-in', str(self.input_file), '-dbtype', self.db_type]
        if self.parse_seqids:
            command.append('-parse_seqids')
        if self.output_db_path != None:
            command +=['-out', str(self.output_db_path)]
        stdout, stderr = run_command(" ".join([str(x) for x in command]))
        print(stdout, stderr)
        return self.output_db_path

class BlastSearch:
    __blast_commands = set(BlastCommands._keys())
    __blast_extensions_nt = frozenset(['.nsq', '.nin', '.nhr'])
    __blast_extensions_pt = frozenset(['.pto', '.ptf', '.phr'])
    __filter_columns = ["qseqid", "sseqid"]

    def __init__(self, db_data: DBData, query_path: Path, blast_params: dict, blast_method: str, blast_columns: List[str], filter_options: Dict[str, FilterOptions]):
        self.db_data = db_data
        if blast_method not in self.__blast_commands:
            raise ValueError("{} is not a valid blast command please pick from {}".format(blast_method, self.__blast_commands))
        self.query_path = query_path
        self.blast_params = blast_params
        self.blast_method = blast_method
        self.blast_columns = blast_columns
        self.filter_options = filter_options

    def get_blast_data(self, db_path: Path, output: Path) -> pd.DataFrame:
        """
        Run blast and parse results
        TODO need to clean up the db_path hand off from the DBData obj its dirty
        """
        
        stdout, stderr = self._run_blast(db=db_path, output=output)
        #logger.info("Blast stdout: {}".format(stdout))
        #logger.info("Blast stderr: {}".format(stderr))
        blast_data = self.parse_blast(output_file=output)
        return blast_data

    def parse_blast(self, output_file: Path):
        """
        Parse a blast output file
        output_file Path: Generate blast data
        """
        df = self.read_hit_table(output_file)
        columns = df.columns.tolist()
        for id_col in self.__filter_columns:
            tp = {}
            if id_col in columns:
                tp[id_col] = 'object'
            df = df.astype(tp)
        for col_name in self.filter_options:
                if col_name in columns:
                    min_value = self.filter_options[col_name].min
                    max_value = self.filter_options[col_name].max
                    include = self.filter_options[col_name].include
                    df = self.filter_df(df,col_name, min_value, max_value, include, columns)
        return df
        

    def filter_df(self,df, col_name,min_value,max_value,include, columns):
        if col_name not in columns:
            return False
        if min_value is not None:
            df = df[df[col_name] >= min_value]
        if max_value is not None:
            df = df[df[col_name] <= max_value]
        if include is not None:
            df = df[df[col_name].isin(include)]
        return df

    def read_hit_table(self, blast_data):
        return pd.read_csv(blast_data,header=None,names=self.blast_columns,sep="\t",low_memory=False)


    def _check_blast_files(self, db_dir: Path, extensions: frozenset):
        """
        """
        extensions_ = set([i.suffix for i in db_dir.iterdir()])
        if not extensions_.issuperset(extensions):
            raise ValueError("Missing required blast files. {}".format([i for i in extensions_ if i not in extensions]))
    
    def validate_blast_db(self, db_data=None):
        """
        """
        if db_data is None:
            db_data = self.db_data
        if db_data.nucleotide:
            self._check_blast_files(db_data.nucleotide, self.__blast_extensions_nt)

        if db_data.protein:
            self._check_blast_files(db_data.protein, self.__blast_extensions_pt)
    
    def _run_blast(self, db: Path, output: Path):
        """
        db PAth: Path to the blast database to use,
        output Path: Path to file for blast output
        """
        command = [
            self.blast_method,
            '-query', self.query_path,
            '-db', str(db),
            '-out', str(output),
            '-outfmt', "'6 {}'".format(' '.join(self.blast_columns)),
        ]
        for param in self.blast_params:
            if param == "parse_seqids":
                command.append(f"-{param}")
            else:
                command += [f'-{param}', f'{self.blast_params[param]}']  
        return run_command(" ".join([str(x) for x in command]))




