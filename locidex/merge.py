import gzip
import json
import os
import re
import sys
import errno
import typing
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
from functools import partial
from mimetypes import guess_type
from multiprocessing import Pool, cpu_count
import logging
import pandas as pd
from locidex.classes.aligner import align, parse_align
from locidex.constants import DBConfig, raise_file_not_found_e
from locidex.report import ReportData, Data, Parameters
from locidex.version import __version__

logger = logging.getLogger(__name__)
logging.basicConfig(filemode=sys.stderr, level=logging.INFO)

class OutputRow(typing.NamedTuple):
    sample_id: str
    keys: list[str]
    message: str

def add_args(parser=None):
    """
    TODO disabling alignment until test data is prepared
    """
    if parser is None:
        parser = ArgumentParser(
            description="Locidex merge: Concatonate set of input profile.json files into  a tsv table or aligned fasta")
    parser.add_argument('-i','--input', type=str, required=True,help='Input file to report', action='append', nargs='+')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file to put results')
    #parser.add_argument('--n_threads','-t', type=int, required=False,
    #                    help='CPU Threads to use', default=1)
    #parser.add_argument('--linker','-l', type=str, required=False,
    #                    help='Linker sequence for alignment', default='NNNNNNNNNNNNNNNNNNNN')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-s', '--strict', required=False, help='Only merge data produces by the same db',
                        action='store_true')
    #parser.add_argument('-a', '--align', required=False, help='Perform alignment with individual loci to produce a concatenated alignment',
    #                    action='store_true')
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    parser.add_argument('-p', '--profile_ref', type=str, required=False, help='Two column TSV file with profile references for overriding MLST profiles. Columns [sample,mlst_alleles]')
    return parser


def get_file_list(input_files):
    file_list = []
    if len(input_files) > 1:
        file_list = input_files
    else:
        if re.search(".json$", input_files[0]) or re.search(".json.gz$", input_files[0]):
            file_list = input_files
        else:
            if not os.path.isfile(input_files[0]):
                logger.critical("File {} does not exists".format(input_files[0]))
                raise_file_not_found_e(input_files[0], logger=logger)

            encoding = guess_type(input_files[0])[1]
            _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
            with _open(input_files[0]) as f:
                for line in f:
                    line = line.rstrip()
                    if not os.path.isfile(line):
                        logger.critical("Could not find file: {}".format(line))
                        raise_file_not_found_e(line, logger)
                    file_list.append(line)
    return file_list

def validate_input_file(data_in: dict, filename: str, db_version: str, db_name: str, perform_db_validation: bool, perform_profile_validation: bool, profile_refs_dict:dict) -> tuple[ReportData, str, str, typing.Optional[OutputRow]]:
    """
    Validate input data for usage verifying db_versions, db_names and MLST profiles are the same
    """
    if perform_profile_validation:
        user_provided_profile = profile_refs_dict[filename][0]
        data_in, mlst_report = validate_and_fix_profiles(data_in, user_provided_profile, filename)
    else:
        mlst_report = None

    try:
        sq_data = ReportData.deseriealize(data_in)
    except KeyError:
        logger.critical("Missing fields in configuration required fields in in reported allele file. Fields required: {}".format(ReportData.fields()))
        raise ValueError("Missing fields in configuration required fields in in reported allele file. Fields required: {}".format(ReportData.fields()))

    else:

        if db_version is not None and sq_data.db_info.db_version != db_version and perform_db_validation:
            logger.critical("You are attempting to merge files that were created using different database versions.")
            raise ValueError("You are attempting to merge files that were created using different database versions.")

        if db_name is not None and sq_data.db_info.db_name != db_name and perform_db_validation:
            logger.critical("You are attempting to merge files that have different names.")
            raise ValueError("You are attempting to merge files that have different names. {} {}".format(sq_data.db_info.db_name, db_name))

    return sq_data, sq_data.db_info.db_version, sq_data.db_info.db_name, mlst_report

def check_files_exist(file_list: list[os.PathLike]) -> None:
    """
    Verify that all files to be analyzed exist
    """
    for file in file_list:
        if not os.path.isfile(file):
            logger.critical("Could not find file: {}".format(file))
            raise_file_not_found_e(file, logger)


def read_file_list(file_list, perform_db_validation=False, perform_profile_validation=False):
    records = {}
    db_version = None
    db_name = None
    check_files_exist(file_list)
    # Before we can perform profile validation of the MLST files we need to check that a proper profile key file was provided
    if perform_profile_validation:
        modified_MLST_files = [["sample", "JSON_key", "error_message"]]
        check_files_exist([perform_profile_validation])
        profile_refs_dict = read_samplesheet(perform_profile_validation)
    else:
        modified_MLST_files = None
        profile_refs_dict = None

    for f in file_list:
        encoding = guess_type(f)[1]
        filename = os.path.basename(f)
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(f) as fh:
            data = json.load(fh)
            sq_data, db_version, db_name, mlst_report = validate_input_file(data,
                                                            filename,
                                                            db_version=db_version,
                                                            db_name=db_name,
                                                            perform_db_validation=perform_db_validation,
                                                            perform_profile_validation=perform_profile_validation,
                                                            profile_refs_dict = profile_refs_dict)
            if perform_profile_validation and mlst_report is not None:
                modified_MLST_files.append(mlst_report)
            sample_name = sq_data.data.sample_name
            if records.get(sq_data.data.sample_name) is None:
                records[sample_name] = sq_data
            else:
                logger.critical("Duplicate sample name detected: {}".format(sq_data.data.sample_name))
                raise ValueError("Attempting to merge allele profiles with the same sample name: {}".format(sq_data.data.sample_name))

    return records, modified_MLST_files

def extract_profiles(records):
    profile = {}
    for id in records:
        for sample_name in records[id].data.profile:
            if profile.get(sample_name) is not None:
                logger.critical("Sample {} already exists and will not be added.")
            profile[sample_name] = records[id].data.profile[sample_name]
    return profile

def extract_seqs(records):
    seqs = {}
    for id in records:
        seqs[id] = records[id].data.seq_data
    return seqs

def write_gene_fastas(seq_data,work_dir):
    d = 0
    files = {}
    for id in seq_data:
        record = seq_data[id]
        for seq_id in record:
            if 'locus_name' not in record[seq_id]:
                continue
            locus_name = record[seq_id]['locus_name']
            if locus_name == '':
                continue
            if 'dna_seq' in record[seq_id]:
                seq = record[seq_id]['dna_seq']
            else:
                seq = record[seq_id]['aa_seq']
            out_file = os.path.join(work_dir, f"{locus_name}.fas")
            if not os.path.isfile(out_file):
                oh = open(out_file,'w')
                files[locus_name] = {'file':out_file}
            else:
                if not locus_name in files:
                    oh = open(out_file,'w')
                    files[locus_name] = {'file':out_file}
                else:
                    oh = open(out_file,'a')
            seq_name = f'{id}'
            oh.write(f'>{seq_name}\n{seq}\n')
            oh.close()
            d+=1
    return files

def read_samplesheet(sample_key_file):
    """
    We will be accepting a file with two columns either ['sample', 'mlst_alleles'] or ['sample_name', 'mlst_alleles']
    """
    try:
        df = pd.read_csv(sample_key_file, sep=",", usecols= lambda column: column in {"sample", "sample_name", "mlst_alleles"}, header=0)
        df["mlst_alleles"] = df["mlst_alleles"].apply(lambda x: os.path.basename(x))
        sampledict = df.set_index("mlst_alleles").T.to_dict('list')


    except:
            logging.critical("File {fname} should be CSV that will contain either ['sample', 'mlst_alleles'] or ['sample_name', 'mlst_alleles']".format(fname = sample_key_file))
            raise FileNotFoundError("Incorrect file format")

    return sampledict

def validate_and_fix_profiles(mlst, sample_id, file_name):
    data_key = "data"
    profile_key = "profile"
    # Extract the profile from the json_data
    profile = mlst.get(data_key, {}).get(profile_key, {})
    # Check for multiple keys in the JSON file and define error message
    keys = sorted(profile.keys())
    original_key = keys[0] if keys else None
    # Define a variable to store the match_status (True or False)
    match_status = sample_id in profile
    # Initialize the error message
    MLST_message = None

    if not keys:
        logger.critical(f"{file_name} is missing the 'profile' section or is completely empty!")
        raise ValueError(f"{file_name} is missing the 'profile' section or is completely empty!")

    elif len(keys) > 1:
        # Check if sample_id matches any key
        if not match_status:
            MLST_message = f"No key in the MLST JSON file ({file_name}) matches the specified sample ID '{sample_id}'. The first key '{original_key}' has been forcefully changed to '{sample_id}' and all other keys have been removed."
            # Retain only the specified sample ID
            mlst[data_key][profile_key] = {sample_id: profile[original_key]}
        else:
            MLST_message = f"MLST JSON file ({file_name}) contains multiple keys: {keys}. The MLST JSON file has been modified to retain only the '{sample_id}' entry"
            # Retain only the specified sample_id in the profile
            mlst[data_key][profile_key] = {sample_id: profile[sample_id]}
    elif not match_status:
        MLST_message = f"{sample_id} ID and JSON key in {file_name} DO NOT MATCH. The '{original_key}' key in {file_name} has been forcefully changed to '{sample_id}': User should manually check input files to ensure correctness."
        # Update the JSON file with the new sample ID
        mlst[data_key][profile_key] = {sample_id: profile[original_key]}
        mlst[data_key]["sample_name"] = sample_id

    # Create a report for all the samples that have their profiles modified in the output profile.tsv
    if MLST_message:
        mlst_report = OutputRow(sample_id, keys, MLST_message)
    else:
        mlst_report = None

    # Write the updated JSON data back to the original file
    return mlst, mlst_report

def run_merge(config):
    analysis_parameters = config

    #Input Parameters
    input_files = config['input'][0]
    outdir = config['outdir']

    ###
    # Commented out as these changes will require test data
    # perform_align = config['align']
    # linker_seq = config['linker']
    # n_threads = config['n_threads']
    ###
    force = config['force']
    validate_db = config['strict']
    profile_refs = config['profile_ref']
    if validate_db is None or validate_db == '':
        validate_db = False
    if profile_refs is None or profile_refs == '':
        profile_refs = False

    run_data = {}
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = analysis_parameters

    if os.path.isdir(outdir) and not force:
        logger.critical(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        raise FileExistsError(errno.EEXIST, os.strerror(errno.EEXIST), str(outdir))

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    #perform merge
    file_list = get_file_list(input_files)
    records, modified_MLST_file_list = read_file_list(file_list,perform_db_validation=validate_db, perform_profile_validation=profile_refs)

    #create profile
    df = pd.DataFrame.from_dict(extract_profiles(records), orient='index')
    df.insert(loc=0, column='sample_id', value=df.index.tolist())
    df.to_csv(os.path.join(outdir,'profile.tsv'),index=False,header=True,sep="\t")

    del(df)
    run_data['result_file'] = os.path.join(outdir,"profile.tsv")

    #Write report of all the MLST files with profile mismatch and how MLST profiles with mismatch were modified
    df = pd.DataFrame(modified_MLST_file_list)
    df.to_csv(f'{outdir}/MLST_error_report.csv', index=False, header=False)


    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    with open(os.path.join(outdir,"run.json"),'w' ) as fh:
        fh.write(json.dumps(run_data, indent=4))


def run(cmd_args=None):
    if cmd_args is None:
        parser = add_args()
        cmd_args = parser.parse_args()
    analysis_parameters = vars(cmd_args)


    config = {}
    for p in analysis_parameters:
        if not p in config:
            config[p] = analysis_parameters[p]

    run_merge(config)


# call main function
if __name__ == '__main__':
    run()


