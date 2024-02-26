import gzip
import json
import os
import re
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
from functools import partial
from mimetypes import guess_type

import pandas as pd

from locidex.version import __version__


def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Locidex merge: Concatonate set of input profile.json files into  a tsv table",
        formatter_class=CustomFormatter)
    parser.add_argument('-i','--input', type=str, required=True,help='Input file to report', action='append', nargs='+')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file to put results')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')

    return parser.parse_args()


def get_file_list(input_files):
    file_list = []
    if len(input_files) > 1:
        file_list = input_files
    else:
        if re.search(".json$", input_files[0]) or re.search(".json.gz$", input_files[0]):
            file_list = input_files
        else:
            encoding = guess_type(input_files[0])[1]
            _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
            with _open(input_files[0]) as f:
                for line in f:
                    file_list.append(line.rstrip())
    return file_list

def read_file_list(file_list):
    records = {}
    for f in file_list:
        if not os.path.isfile(f):
            continue
        encoding = guess_type(f)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(f) as fh:
            data = json.load(fh)
            records = records | data
    return records



def run():
    cmd_args = parse_args()
    analysis_parameters = vars(cmd_args)

    #Input Parameters
    input_files = cmd_args.input[0]
    outdir = cmd_args.outdir
    force = cmd_args.force



    run_data = {}
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = vars(cmd_args)

    if os.path.isdir(outdir) and not force:
        print(f'Error {outdir} exists, if you would like to overwrite, then specify --force')
        sys.exit()

    if not os.path.isdir(outdir):
        os.makedirs(outdir, 0o755)

    file_list = get_file_list(input_files)
    records = read_file_list(file_list)

    df = pd.DataFrame.from_dict(records, orient='index')
    df.insert(loc=0,
              column='sample_id',
              value=df.index.tolist())
    df.to_csv(os.path.join(outdir,'profile.tsv'),index=False,header=True,sep="\t")


# call main function
if __name__ == '__main__':
    run()


