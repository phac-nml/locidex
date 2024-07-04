#!/usr/bin/env python3

import sys
import argparse
import traceback
import logging
from . import format, extract, report, merge, search, build, manifest, utils, constants

logger = logging.getLogger(__name__)
logging.basicConfig(filemode=sys.stderr, level=logging.DEBUG)

tasks = {
    'search': (search, 'Query set of Loci/Genes against a database to produce a sequence store for downstream processing'),
    'extract': (extract, 'Extract loci from a genome based on a locidex database'),
    'report': (report, 'Filter a sequence store and produce and extract of blast results and gene profile'),
    'merge':  (merge, 'Merge a set of gene profiles into a standard profile format'),
    'format': (format, 'Format fasta files from other MLST databases for use with locidex build'),
    'build': (build, 'Build a locidex database'),
    'manifest': (manifest, 'Create a multi-database folder manifest'),
}


def main(argv=None):
    module_idx = 0
    help_msg = 1
    parser = argparse.ArgumentParser(prog="locidex")
    sub_parsers = parser.add_subparsers(dest="command")
    for k, v in tasks.items():
        format_parser = sub_parsers.add_parser(k, description=v[help_msg], help=v[help_msg])
        v[module_idx].add_args(format_parser)
    
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        try:
            utils.check_utilities(logger, constants.UTILITIES_CHECK)
        except FileNotFoundError:
            pass
        sys.exit()

    error_file = "errors.txt"
    try:
        utils.check_utilities(logger, constants.UTILITIES_CHECK)
        logger.info("Running {}".format(args.command))
        tasks[args.command][module_idx].run(args)
        logger.info("Finished: {}".format(args.command))
    except Exception as e:
        with open(error_file, "w") as f:
            f.write(traceback.format_exc())
        error_number = e.errno if hasattr(e, "errno") else -1
        logger.critical("Program exited with errors, please review logs. For the full traceback please see file: {}".format(error_file))
        SystemExit(error_number)
    else:
        logger.info("Program finished without errors.")

# call main function
if __name__ == '__main__':
    main()