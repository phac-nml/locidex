import pathlib
import json
import os
import re
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from datetime import datetime
from locidex.version import __version__


def add_args(parser=None):
    if parser is None:
        parser = ArgumentParser(
            description="Locidex manifest: Setup directory of databases for use with search")
    parser.add_argument('-i','--input', type=str, required=True,help='Input directory of locidex databases')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)
    return parser

def run_merge(config):
    analysis_parameters = config

    #Input Parameters
    input_dir = config['input']
    in_dirname = input_dir.split('/')[-1]

    run_data = {}
    run_data['analysis_start_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    run_data['parameters'] = analysis_parameters

    db_keys = [
        "db_name",
        "db_version",
        "db_date",
        "db_author",
        "db_desc",
        "db_num_seqs",
    ]

    d = pathlib.Path(input_dir).rglob('*')
    config_files = {}
    for item in d:
        if item.is_dir():
            continue
        fpath = item.resolve()
        dirname = os.path.dirname(fpath).split('/')[-1]
        fname = os.path.basename(item)
        if fname != 'config.json':
            continue
        c = {}
        with open(fpath ,'r') as fh:
            c = json.load(fh)
        if len(c) == 0:
            continue
        for field in db_keys:
            if not field in c:
                print(f'Error db config: {fpath} is missing a needed field key for {field}, please set one', file=sys.stderr)
                raise KeyError

            v = c[field]
            if v == '':
                print(f'Error db config: {fpath} is missing a needed field value for {field}, please set one', file=sys.stderr)
                raise KeyError
        
        db_name = str(c['db_name'])
        db_version = str(c['db_version'])
        if not db_name in config_files:
            config_files[db_name] = {}
        if db_version in config_files[db_name]:
            print(f"Error you are trying to populate duplicate entries for db_name {db_name} and version {db_version}. \
                  Manifest only supports distinct db_entries, please resolve duplicates", file=sys.stderr)
            sys.exit()
        
        config_files[db_name][db_version] = {
            #'db_relative_path_dir': f"{in_dirname}/{dirname}",
            'db_relative_path_dir': os.path.join(in_dirname, dirname),
            #'db_relative_path_config': f"{in_dirname}/{dirname}/config.json",
            'db_relative_path_config': os.path.join(in_dirname, dirname, "config.json"),
        }

    with open(os.path.join(input_dir,"manifest.json"),'w' ) as fh:
        fh.write(json.dumps(config_files, indent=4))

    run_data['analysis_end_time'] = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    with open(os.path.join(input_dir,"run.json"),'w' ) as fh:
        fh.write(json.dumps(run_data, indent=4))



def run(cmd_args=None):
    #cmd_args = parse_args()
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


