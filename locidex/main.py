#!/usr/bin/env python3

import sys
import argparse
from . import format, extract, report, merge, search, build, manifest

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
        sys.exit()
    #print("args", tasks[args.command][module_idx].run(args))
    tasks[args.command][module_idx].run(args)


# call main function
if __name__ == '__main__':
    sys.exit(main())