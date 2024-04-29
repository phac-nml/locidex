import pathlib
import json
from typing import List, Union, Tuple, Dict
import os
import re
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
from locidex.version import __version__
from locidex.constants import DBConfig, DBFiles, ManifestFields


def add_args(parser=None):
    if parser is None:
        parser = ArgumentParser(
            description="Locidex manifest: Setup directory of databases for use with search")
    parser.add_argument('-i','--input', type=str, required=True,help='Input directory containing multiplie locidex databases')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    return parser


def check_config(directory: pathlib.Path) -> None:
    """
    Validate config file in a directory. Throws an error if any required parameters
    are missing.

    directory: Path of the directory containing the parent.
    """

    config_dir = pathlib.Path(directory / DBFiles.config_file)
    config_data: Union[DBConfig, None] = None 
    with open(config_dir, 'r') as conf:
        config_data = DBConfig(**json.load(conf))
        for k, v in config_data.to_dict().items():
            if v is None or v == '':
                raise AttributeError("Config cannot have missing values: {}".format(k))
    return config_data

def validate_db_files(allele_dir: List[pathlib.Path]) -> List[Tuple[pathlib.Path, DBConfig]]:
    """
    Validates a directory of allele databases, and verifies that the config contains
    the required fields

    allele_dir List[pathlib.Path: Directory of various allele databases needed by mikrokondo
    """
    db_configs: Tuple[pathlib.Path, DBConfig] = []
    for a_dir in allele_dir:
        for k, v in DBFiles.items():
            if not pathlib.Path(a_dir / v).exists():
                raise FileNotFoundError("Required file {} does not exist.".format(k))
        db_configs.append((a_dir, check_config(a_dir)))
    return db_configs


def check_dbs(file_in: pathlib.Path) -> List[pathlib.PosixPath]:
    """
    Checks that all locidex databases in a directory are complete.

    file_in: A path to a directory of databases
    """
    db_dirs = [p for p in file_in.iterdir() if p.is_dir()]
    return db_dirs

def create_manifest(file_in: pathlib.Path):
    """
    Create a manifest file for each of the locidex dbs.

    file_in pathlib.Path: File path to directory of databases
    """
    allele_dirs: List[pathlib.Path] = check_dbs(file_in)
    validated_dbs: List[Tuple[pathlib.Path, DBConfig]] = validate_db_files(allele_dirs)
    db_manifest = dict()
    for path, conf in validated_dbs:
        if db_manifest.get(conf.db_name) is not None:
            raise KeyError("Databases with the same name have been specified (name: {}, path: {})".format(conf.db_name, path))

        db_manifest[conf.db_name] = {
            ManifestFields.db_path: str(path),
            ManifestFields.config_data: conf.to_dict() 
        }
    return db_manifest


#def run_merge(config):
#    analysis_parameters = config
#
#    #Input Parameters
#    input_dir = config['input']
#    in_dirname = input_dir.split('/')[-1]
#
#    run_data = {}
#    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
#    run_data['parameters'] = analysis_parameters
#
#    db_keys = DBConfig.keys()
#
#    d = pathlib.Path(input_dir).rglob('*')
#    config_files = {}
#    for item in d:
#        if item.is_dir():
#            continue
#        fpath = item.resolve()
#        dirname = os.path.dirname(fpath).split('/')[-1]
#        fname = os.path.basename(item)
#        if fname != 'config.json':
#            continue
#        c = {}
#        with open(fpath ,'r') as fh:
#            c = json.load(fh)
#        if len(c) == 0:
#            continue
#        for field in db_keys:
#            if not field in c:
#                print(f'Error db config: {fpath} is missing a needed field key for {field}, please set one', file=sys.stderr)
#                raise KeyError
#
#            v = c[field]
#            if v == '':
#                print(f'Error db config: {fpath} is missing a needed field value for {field}, please set one', file=sys.stderr)
#                raise KeyError
#        
#        db_name = str(c['db_name'])
#        db_version = str(c['db_version'])
#        if not db_name in config_files:
#            config_files[db_name] = {}
#        if db_version in config_files[db_name]:
#            print(f"Error you are trying to populate duplicate entries for db_name {db_name} and version {db_version}. \
#                  Manifest only supports distinct db_entries, please resolve duplicates", file=sys.stderr)
#            sys.exit()
#        
#        config_files[db_name][db_version] = {
#            #'db_relative_path_dir': f"{in_dirname}/{dirname}",
#            'db_relative_path_dir': os.path.join(in_dirname, dirname),
#            #'db_relative_path_config': f"{in_dirname}/{dirname}/config.json",
#            'db_relative_path_config': os.path.join(in_dirname, dirname, "config.json"),
#        }
#
#    with open(os.path.join(input_dir,"manifest.json"),'w' ) as fh:
#        fh.write(json.dumps(config_files, indent=4))
#
#    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
#    with open(os.path.join(input_dir,"run.json"),'w' ) as fh:
#        fh.write(json.dumps(run_data, indent=4))
#
#
#
#def run(cmd_args=None):
#    #cmd_args = parse_args()
#    if cmd_args is None:
#        parser = add_args()
#        cmd_args = parser.parse_args()
#    analysis_parameters = vars(cmd_args)
#
#
#    config = {}
#    for p in analysis_parameters:
#        if not p in config:
#            config[p] = analysis_parameters[p]
#
#    run_merge(config)
#
#
## call main function
#if __name__ == '__main__':
#    run()


