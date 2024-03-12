#!/usr/bin/env python3

import sys

tasks = {
    'search': 'Query set of Loci/Genes against a database to produce a sequence store for downstream processing',
    'extract': 'Extract loci from a genome based on a locidex database',
    'report': 'Filter a sequence store and produce and extract of blast results and gene profile',
    'merge':  'Merge a set of gene profiles into a standard profile format',
    'format': 'Format fasta files from other MLST databases for use with locidex build',
    'build': 'Builds locidex db',
}



ordered_tasks = [
    'search',
    'extract',
    'report',
    'merge',
    'format',
    'build'
]


def print_usage_and_exit():
    """
    This method prints brief usage instructions of Clade-o-matic to the command line
    """
    print('Usage: locidex <command> [options] <required arguments>', file=sys.stderr)
    print('\nTo get minimal usage for a command use:\ngas command', file=sys.stderr)
    print('\nTo get full help for a command use one of:\nlocidex command -h\nlocidex command --help\n', file=sys.stderr)
    print('\nAvailable commands:\n', file=sys.stderr)
    max_task_length = max([len(x) for x in list(tasks.keys())]) + 1
    for task in ordered_tasks:
        print('{{0: <{}}}'.format(max_task_length).format(task), tasks[task], sep=' ', file=sys.stderr)
    sys.exit(0)

def main():

    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '-help', '--help']:
        print_usage_and_exit()

    task = sys.argv.pop(1)

    if task not in tasks:
        print('Task "' + task + '" not recognised. Cannot continue.\n', file=sys.stderr)
        print_usage_and_exit()

    exec('import locidex.' + task)
    exec('locidex.' + task + '.run()')

# call main function
if __name__ == '__main__':
    main()