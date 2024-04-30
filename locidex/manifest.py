import pathlib
import json
from typing import List, Union, Tuple, Dict
from dataclasses import dataclass
import os
import re
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
from locidex.version import __version__
from locidex.constants import DBConfig, DBFiles




class ManifestItem:
    """
    Manifest item created for exporting and importing locidex items
    """
    __path_key = 'path'
    __config_key = 'config'

    def __init__(self, db: pathlib.Path, root_db: pathlib.Path, config: DBConfig):
        """
        db Path: Relative path to the allele database
        root_db Path: The directory containing the manifest.json file required to resolve paths
        config DBConfig: Database configuration data
        """
        self.db = db
        self.config = config
        self.root_db = root_db

    @property
    def db_path(self):
        return self.root_db / self.db

    def to_dict(self):
        return {self.__path_key: str(self.db), self.__config_key: self.config.to_dict()}
    
    def __repr__(self) -> str:
        return "Allele location: {}\n Config data: {}".format(self.db, self.config)
    
    @classmethod
    def path_key(cls):
        """
        ! Not passing this as a property for the class method as apparently that will be deprecated in 3.13
        """
        return cls.__path_key
    
    @classmethod
    def config_key(cls):
        """
        ! Not passing this as a property for the class method as apparently that will be deprecated in 3.13
        """
        return cls.__config_key


@dataclass(frozen=True)
class _Constants:
    manifest_name: pathlib.Path = "manifest.json"

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

def validate_db_files(allele_dir: List[pathlib.Path], file_in: pathlib.Path) -> List[Tuple[pathlib.Path, DBConfig]]:
    """
    Validates a directory of allele databases, and verifies that the config contains
    the required fields

    allele_dir List[pathlib.Path: Directory of various allele databases needed by locidex
    file_in Path: Root directory to set files relative too
    """
    db_configs: Tuple[pathlib.Path, DBConfig] = []
    for a_dir in allele_dir:
        for k, v in DBFiles.items():
            if not pathlib.Path(a_dir / v).exists():
                raise FileNotFoundError("Required file {} does not exist.".format(k))
        db_configs.append((a_dir.relative_to(file_in), check_config(a_dir)))
    return db_configs


def check_dbs(file_in: pathlib.Path) -> List[pathlib.PosixPath]:
    """
    Checks that all locidex databases in a directory are complete.

    file_in: A path to a directory of databases
    """
    db_dirs = [p for p in file_in.iterdir() if p.is_dir()]
    return db_dirs

def create_manifest(file_in: pathlib.Path) -> Dict[str, List[Dict[str, str]]]:
    """
    Create a manifest file for each of the locidex dbs.

    file_in pathlib.Path: File path to directory of databases
    """
    allele_dirs: List[pathlib.Path] = check_dbs(file_in)
    validated_dbs: List[Tuple[pathlib.Path, DBConfig]] = validate_db_files(allele_dirs, file_in)
    db_manifest = dict()
    for path, conf in validated_dbs:

        if db_manifest.get(conf.db_name) is None:
            db_manifest[conf.db_name] = []
            
        if db_manifest[conf.db_name] and (versions := [i.db_version for i in db_manifest[conf.db_name]]):
            if conf.db_version in versions:
                raise KeyError("Databases with the same name and version have been specified (name: {}, path: {}, version: {})".format(conf.db_name, path, conf.db_version))
            
        db_manifest[conf.db_name].append(ManifestItem(db=path, config=conf, root_db=file_in).to_dict())
    return db_manifest


def write_manifest(file_in: pathlib.Path, manifest: Dict[str, List[Dict[str, str]]]) -> pathlib.Path:
    """
    Write the manifest.json file

    file_in Path: Specified input directory
    manifest dict: data to write to the manifest
    """

    manifest_file = _Constants.manifest_name
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

def select_db(manifest_data: Dict[str, List[ManifestItem]], name: str, version: str):
    """
    Select a locidex database from the manifest file provided.

    manifest_data Dict[str, List[ManifestItem]]: Parsed manifest file data for selecting a database
    name str: Name of database to select
    version str: version of selected database to select
    """
    db_data = manifest_data.get(name)
    if db_data is None:
        raise KeyError("Could not find database with specified name: {}".format(name))
    
    try:
        db = next(filter(lambda x: x.config.db_version == version, db_data))
    except StopIteration:
        raise ValueError("No database entry with version: {}".format(version))
    
    return db

def read_manifest(input_file: pathlib.Path) -> dict:
    """
    input_file Path: Manifest file to be parsed
    """
    if not input_file.is_dir():
        raise AssertionError("Allele database directory must be passed directly.")
    
    manifest_file = input_file / _Constants.manifest_name
    manifest_data: Dict[str, ManifestItem] = dict()
    with open(manifest_file, 'r', encoding='utf8') as mani_in:
        manifest = json.load(mani_in)
        for k, list_manifests in manifest.items():
            if manifest_data.get(k) is None:
                manifest_data[k] = []
            for v in list_manifests:
                manifest_item = ManifestItem(db=v[ManifestItem.path_key()], config=DBConfig(**v[ManifestItem.config_key()]), root_db=input_file)
                manifest_data[k].append(manifest_item)
    return manifest_data

def get_manifest_db(input_file: pathlib.Path, name: str, version: str):
    output = read_manifest(input_file)
    db_out = select_db(output, name, version)
    return db_out.db_path

# call main function
if __name__ == '__main__':
    run()

