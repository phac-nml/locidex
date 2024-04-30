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


def write_manifest(file_in: pathlib.Path, manifest: Dict[str, Dict[str, Union[str, Dict[str, str]]]]) -> None:
    """
    Write the manifest.json file

    file_in Path: Specified input directory
    manifest dict: data to write to the manifest
    """

    manifest_file = "manifest.json"
    path_out = file_in.joinpath(manifest_file)
    with open(path_out, 'w', encoding='utf8') as m_out:
        json.dump(manifest, m_out, indent=2)
    return path_out


def run(cmd_args):
    if cmd_args is None:
        parser = add_args()
        cmd_args = parser.parse_args()
    directory_in = pathlib.Path(cmd_args.input)
    directory_in.exists()
    manifest = create_manifest(directory_in)
    return write_manifest(directory_in, manifest)


# call main function
if __name__ == '__main__':
    run()

