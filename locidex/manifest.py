import pathlib
import json
from typing import List, Union, Tuple, Dict, Optional
from dataclasses import dataclass
import os
import re
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
from locidex.version import __version__
from locidex.constants import DBConfig, DBFiles, raise_file_not_found_e
import logging
import errno

logger = logging.getLogger(__name__)
logging.basicConfig(filemode=sys.stderr, level=logging.DEBUG)

class DBData:
    """
    Validate and get all database data for other modules.
    
    * This class will create some redundancy and will need to be refactored to reflect
    * the overall use of this module in the future. But once refactoring is complete 
    * it should be much easier to refactor this and the other modules. Additionally, at that
    * point we should have a better understanding of how all the modules fit together.
    """

    __nucleotide_name = "nucleotide"
    __nucleotide_db_type = "nucl"
    __protein_name = "protein"
    __protein_db_type = "prot"
    __nucleotide_path = pathlib.Path(__nucleotide_name)
    __protein_path = pathlib.Path(__protein_name)

    def __init__(self, db_dir: pathlib.Path):
        self.db_dir = pathlib.Path(db_dir)
        self.config_data: DBConfig = self._get_config(self.db_dir)
        self.metadata: dict = self._get_metadata(self.db_dir)
        self.nucleotide, self.protein = self._get_blast_dbs(self.db_dir, self.config_data)

    @classmethod
    def nucleotide_db_type(cls):
        return cls.__nucleotide_db_type

    @classmethod
    def protein_db_type(cls):
        return cls.__protein_db_type

    @classmethod
    def protein_name(cls):
        return cls.__protein_name
    
    @classmethod
    def nucleotide_name(cls):
        return cls.__nucleotide_name

    @property
    def nucleotide_blast_db(self):
        if self.nucleotide is None:
            logger.critical("Nucleotide blast database does not exist.")
            raise_file_not_found_e(self.nucleotide, logger)
        return self.nucleotide / self.__nucleotide_path
    
    @property
    def protein_blast_db(self):
        if self.protein is None:
            logger.critical("Protein blast database does not exist.")
            raise_file_not_found_e(self.protein, logger)
        return self.protein / self.__protein_path

    def _get_config(self, db_dir: pathlib.Path) -> DBConfig:
        """
        Validates the config file and searializes the data into a DBConfig object
        """
        return check_config(db_dir)
    
    def _get_metadata(self, db_dir: pathlib.Path) -> dict:
        metadata_file = db_dir.joinpath(DBFiles.meta_file)
        if not metadata_file.exists():
            logger.critical("Metadata file does not appear to exist in db: {}".format(db_dir))
            raise_file_not_found_e(str(metadata_file), logger)
        md_data = None
        with open(metadata_file, 'r') as md:
            md_data = json.load(md)
        return md_data


    def _get_blast_dbs(self, db_dir: pathlib.Path, config_data: DBConfig) -> Tuple[Optional[pathlib.Path], Optional[pathlib.Path]]:
        blast_db = db_dir.joinpath(DBFiles.blast_dir)
        nucleotide: Optional[pathlib.Path] = None
        protein: Optional[pathlib.Path] = None
        if not blast_db.exists():
            logger.critical("blast directory not found. Database path maybe incorrect: {}".format(str(db_dir)))
            raise NotADirectoryError(errno.ENOTDIR, os.strerror(errno.ENOTDIR), str(db_dir))
        if config_data.is_nucl:
            nucleotide = blast_db.joinpath(self.__nucleotide_path)
            if not nucleotide.exists():
                logger.critical("Cannot find nucleotide database, but it should exist. {}".format(nucleotide))
                raise_file_not_found_e(nucleotide, logger)
        if config_data.is_prot:
            protein = blast_db.joinpath(self.__protein_path)
            if not protein.exists():
                logger.critical("Cannot find protein database, but it should exist. {}".format(protein))
                raise_file_not_found_e(protein, logger)
        return nucleotide, protein


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
    parser.add_argument('-i','--input', type=str, required=True,help='Input directory containing multiple locidex databases')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    return parser


def check_config(directory: pathlib.Path) -> DBConfig:
    """
    Validate config file in a directory. Throws an error if any required parameters
    are missing.

    directory: Path of the directory containing the parent.
    """

    config_dir = pathlib.Path(directory).joinpath(DBFiles.config_file)
    config_data: Optional[DBConfig] = None
    if not config_dir.exists():
        logger.critical("Could not find config file: {}".format(config_dir))
        raise_file_not_found_e(config_dir, logger)
    
    with open(config_dir, 'r') as conf:
        config_data = DBConfig(**json.load(conf))
        for k, v in config_data.to_dict().items():
            if v is None or v == '':
                logger.critical("Missing value in config file for key {} which has value {}".format(k, v))
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
                logger.critical("Required file {} does not exist.".format(k))
                raise_file_not_found_e(str(a_dir / v), logger)
            
        db_configs.append((a_dir.relative_to(file_in), check_config(a_dir)))
    return db_configs


def check_dbs(file_in: pathlib.Path) -> List[pathlib.PosixPath]:
    """
    Checks that all locidex databases in a directory are complete.

    file_in: A path to a directory of databases
    """
    logger.debug("Checking that the following databases exist in: {}".format(file_in))
    db_dirs = [p for p in file_in.iterdir() if p.is_dir()]
    if not db_dirs:
        logger.critical("No valid databases found in: {}".format(file_in))
        raise AssertionError("No valid databases found in: {}".format(file_in))
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
        logger.info("Adding database: {} to manifest.".format(str(path)))
        if db_manifest.get(conf.db_name) is None:
            db_manifest[conf.db_name] = []
            
        if db_manifest[conf.db_name] and (versions := [i.db_version for i in db_manifest[conf.db_name]]):
            if conf.db_version in versions:
                logger.critical("Databases with the same name and version have been specified (name: {}, path: {}, version: {})".format(conf.db_name, path, conf.db_version))
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
    logger.debug("Creating manifest file in directory: {} with name: {}".format(file_in, manifest_file))
    path_out = file_in.joinpath(manifest_file)
    logger.info("Writing manifest file to: {}".format(str(path_out)))
    with open(path_out, 'w', encoding='utf8') as m_out:
        json.dump(manifest, m_out, indent=2)
    return path_out


def run(cmd_args=None):
    if cmd_args is None:
        parser = add_args()
        cmd_args = parser.parse_args()
    directory_in = pathlib.Path(cmd_args.input)
    if not directory_in.exists():
        logger.critical("Directory: {} does not appear to exist.".format(str(directory_in)))
        raise NotADirectoryError(errno.ENOTDIR, os.strerror(errno.ENOTDIR), str(directory_in))

    manifest = create_manifest(directory_in)
    return write_manifest(directory_in, manifest)

def select_db(manifest_data: Dict[str, List[ManifestItem]], name: str, version: str) -> ManifestItem:
    """
    Select a locidex database from the manifest file provided.

    manifest_data Dict[str, List[ManifestItem]]: Parsed manifest file data for selecting a database
    name str: Name of database to select
    version str: version of selected database to select
    """
    db_data = manifest_data.get(name)
    if db_data is None:
        logger.critical("Could not find database with specified name: {}".format(name))
        raise KeyError("Could not find database with specified name: {}".format(name))
    
    try:
        db = next(filter(lambda x: x.config.db_version == version, db_data))
    except StopIteration:
        logger.critical("No database entry with version: {}".format(version))
        raise ValueError("No database entry with version: {}".format(version))
    
    return db

def read_manifest(input_file: pathlib.Path) -> dict:
    """
    input_file Path: Manifest file to be parsed
    """
    if not input_file.is_dir():
        logger.critical("Please pass the database directory, not a file.")
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

def get_manifest_db(input_file: pathlib.Path, name: str, version: str) -> ManifestItem:
    """
    Returns path to the database file selected
    """
    output = read_manifest(input_file)
    db_out = select_db(output, name, version)
    return db_out

# call main function
if __name__ == '__main__':
    run()

