import json
import os
import sys
from argparse import (
    ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
from mimetypes import guess_type
from functools import partial
import gzip
import logging
from locidex.version import __version__

logger = logging.getLogger(__name__)


def add_args(parser=None):
    if parser is None:
        parser = ArgumentParser(
            description="Locidex: Advanced searching and filtering of sequence databases using query sequences",)
    parser.add_argument('-p', '--profile', type=str,
                        required=True, help='Integer formatted profile')
    parser.add_argument(
        "--mapping_file",
        "-m",
        type=str,
        required=True,
        help="json formatted allele mapping",
    )
    parser.add_argument('-o', '--outfile', type=str,
                        required=True, help='Output profile to put results')
    parser.add_argument('-V', '--version', action='version',
                        version="%(prog)s " + __version__)
    parser.add_argument('-f', '--force', required=False, help='Overwrite existing directory',
                        action='store_true')
    return parser


def parse_allele_map(allele_mapping_file,reverse=False):
    if not os.path.isfile(allele_mapping_file):
        return {}
    encoding = guess_type(all)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with _open(allele_mapping_file, "r", encoding="utf-8") as mapping_fh:
        allele_map = json.loads(mapping_fh.read())

    if len(allele_map) == 0 or reverse:
        return allele_map

    locus_lookup = {}
    loci = list(allele_map.keys)
    for l in loci: 
        locus_lookup[l] = {}
        for hdx in allele_map[l]:
            idx = str(allele_map[l][hdx])
            locus_lookup[l][idx] = hdx
        del (allele_map[l])
    
    return locus_lookup

def list_to_string(data, sep="\t" ):
    return f'{sep}'.join([str(x) for x in data])

def hashify_profile(profile_file, out_file, locus_lookup, sep="\t",file_encoding='utf-8'):
    encoding = guess_type(all)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with open(out_file,'w') as oh:
        with _open(profile_file, "r", encoding=file_encoding) as fh:
            header = fh.read().rstrip().split("\t")
            oh.write(f'{list_to_string(header,sep=sep)}\n')
            for line in fh:
                line = [str(x) for x in line.split(sep)]
                for idx, locus_id in enumerate(header):
                    allele_num = line[idx]
                    if locus_id in locus_lookup:
                        if allele_num in locus_lookup[locus_id]:
                            line[idx] = locus_lookup[locus_id][allele_num]
                oh.write(f'{list_to_string(line,sep=sep)}\n')

def run_profile_conversion(hash=True):
    return

def run(cmd_args=None):
    logger.info("Beginning conversion from integer to hash based profiles")
    if cmd_args is None:
        parser = add_args()
        cmd_args = parser.parse_args()
    analysis_parameters = vars(cmd_args)
    profile_file = cmd_args.profile

    run_profile_conversion(profile_file)
    logger.info("Finishing report workflow.")


# call main function
if __name__ == '__main__':
    run()








