import gzip
import json
import os
import re
import sys
import errno
import csv
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
    parser.add_argument('-k', '--samplekey', type=str, required=False, help='Two column TSV file for overriding MLST profile key of JSON file. Columns [sample,mlst_alleles]')
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

def validate_input_file(data_in: dict, db_version: str, db_name: str, perform_validation: bool) -> tuple[ReportData, str, str]:
    """
    Validate input data for usage verifying db_versions and db_names are the same
    """

    try:
        sq_data = ReportData.deseriealize(data_in)
    except KeyError:
        logger.critical("Missing fields in configuration required fields in in reported allele file. Fields required: {}".format(ReportData.fields()))
        raise ValueError("Missing fields in configuration required fields in in reported allele file. Fields required: {}".format(ReportData.fields()))

    else:

        if db_version is not None and sq_data.db_info.db_version != db_version and perform_validation:
            logger.critical("You are attempting to merge files that were created using different database versions.")
            raise ValueError("You are attempting to merge files that were created using different database versions.")

        if db_name is not None and sq_data.db_info.db_name != db_name and perform_validation:
            logger.critical("You are attempting to merge files that have different names.")
            raise ValueError("You are attempting to merge files that have different names. {} {}".format(sq_data.db_info.db_name, db_name))

    return sq_data, sq_data.db_info.db_version, sq_data.db_info.db_name

def check_files_exist(file_list: list[os.PathLike]) -> None:
    """
    Verify that all files to be analyzed exist
    """
    for file in file_list:
        if not os.path.isfile(file):
            logger.critical("Could not find file: {}".format(file))
            raise_file_not_found_e(file, logger)


def read_file_list(file_list, outputdirectory, perform_validation=False, key_sample_name=None):
    records = {}
    db_version = None
    db_name = None
    check_files_exist(file_list)
    error_reports = []

    for f in file_list:
        if key_sample_name:
            alt_profile = key_sample_name[os.path.basename(f)][0]
        encoding = guess_type(f)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(f) as fh:
            data = json.load(fh)
            if key_sample_name:
                data, compare_errmsg = compare_profiles(data,alt_profile, os.path.basename(f))
                if compare_errmsg:
                    error_reports.append(compare_errmsg)
                        # Write the a new updated JSON data back to a new file
                    with gzip.open("{}/{}.gz".format(outputdirectory, os.path.basename(f)), "wt") as f:
                        json.dump(data, f, indent=4)
            sq_data, db_version, db_name = validate_input_file(data,
                                                            db_version=db_version,
                                                            db_name=db_name,
                                                            perform_validation=perform_validation)

            sample_name = sq_data.data.sample_name
            if records.get(sq_data.data.sample_name) is None:
                records[sample_name] = sq_data
            else:
                logger.critical("Duplicate sample name detected: {}".format(sq_data.data.sample_name))
                raise ValueError("Attempting to merge allele profiles with the same sample name: {}".format(sq_data.data.sample_name))

    return records, error_reports

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
        f =open(sample_key_file,'r')
        sampledict = {}
        for line in f:
            line = line.rstrip().split(",")
            if len(line) != 2:
                logger.critical("File should be tsv with two columns [sample,mlst_alleles]")
                raise_file_not_found_e(logger)
            sample = line[0]
            mlst_file = os.path.basename(line[1])
            sampledict[mlst_file] = [sample]
        return sampledict

def compare_profiles(mlst, sample_id, file_name):
    # Extract the profile from the json_data
    profile = mlst.get("data", {}).get("profile", {})
    # Check for multiple keys in the JSON file and define error message
    keys = sorted(profile.keys())
    original_key = keys[0] if keys else None
    # Define a variable to store the match_status (True or False)
    match_status = sample_id in profile
    # Initialize the error message
    error_message = None

    if not keys:
        logger.critical(f"{file_name} is missing the 'profile' section or is completely empty!")
        raise ValueError(f"{file_name} is missing the 'profile' section or is completely empty!")

    elif len(keys) > 1:
        # Check if sample_id matches any key
        if not match_status:
            error_message = f"No key in the MLST JSON file ({file_name}) matches the specified sample ID '{sample_id}'. The first key '{original_key}' has been forcefully changed to '{sample_id}' and all other keys have been removed."
            # Retain only the specified sample ID
            mlst["data"]["profile"] = {sample_id: profile.pop(original_key)}
        else:
            error_message = f"MLST JSON file ({file_name}) contains multiple keys: {keys}. The MLST JSON file has been modified to retain only the '{sample_id}' entry"
            # Retain only the specified sample_id in the profile
            mlst["data"]["profile"] = {sample_id: profile[sample_id]}
    elif not match_status:
        error_message = f"{sample_id} ID and JSON key in {file_name} DO NOT MATCH. The '{original_key}' key in {file_name} has been forcefully changed to '{sample_id}': User should manually check input files to ensure correctness."
        # Update the JSON file with the new sample ID
        mlst["data"]["profile"] = {sample_id: profile.pop(original_key)}
        mlst["data"]["sample_name"] = sample_id

    error_report = [sample_id, keys, error_message]
    # Write the updated JSON data back to the original file
    return mlst, error_report

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
    samplekeys = config['samplekey']
    if validate_db is None or validate_db == '':
        validate_db = False

    run_data = {}
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = analysis_parameters

    if os.path.isdir(outdir) and not force:
        logger.critical(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        raise FileExistsError(errno.EEXIST, os.strerror(errno.EEXIST), str(outdir))

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    if samplekeys:
        check_files_exist([samplekeys])
        sample_dict = read_samplesheet(samplekeys)
    else:
        sample_dict = None

    #perform merge
    file_list = get_file_list(input_files)
    records, compare_error = read_file_list(file_list, outdir, perform_validation=validate_db, key_sample_name=sample_dict)

    #create profile
    df = pd.DataFrame.from_dict(extract_profiles(records), orient='index')
    df.insert(loc=0, column='sample_id', value=df.index.tolist())
    df.to_csv(os.path.join(outdir,'profile.tsv'),index=False,header=True,sep="\t")

    del(df)
    run_data['result_file'] = os.path.join(outdir,"profile.tsv")

    #Write error messages for profile mismatch (compare_profiles())
    for error_message in compare_error:
            if error_message[2]:
                output_error_file = outdir + "/" + error_message[0] + "_error_report.csv"
                with open(output_error_file, "w", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["sample", "JSON_key", "error_message"])
                    writer.writerow([error_message[0], error_message[1], error_message[2]])

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


