#!/usr/bin/env python3
import os
from distutils.core import setup

from setuptools import find_packages

from locidex.version import __version__

author = 'James Robertson'

classifiers = """
Development Status :: 4 - Beta
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


exec(open('locidex/version.py').read())

setup(
    name='locidex',
    include_package_data=True,
    version=__version__,
    python_requires='>=3.8.2,<4',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    packages=find_packages(exclude=['tests']),
    url='https://github.com/phac-nml/locidex',
    license='GPLv3',
    author='James Robertson',
    author_email='james.robertson@phac-aspc.gc.ca',
    description=(
        'Genomic Address Service: De novo clustering and cluster address assignment'),
    keywords='cgMLST, wgMLST, outbreak, surveillance, blast, antimicrobial resistance, amr, virulence',
    classifiers=classifiers,
    package_dir={'locidex': 'locidex'},
    package_data={
        "": ["*.txt","*.fas","*.fasta","*.json","nucleotide*","protein*","*.gbk"],
    },

    install_requires=[
        'numba>=0.57.1',
        'numpy>=1.24.4',
        'tables>=3.8.0',
        'six>=1.16.0',
        'pandas>=2.0.2 ',
        'biopython>=1.83',
        'pyrodigal'

    ],

    entry_points={
        'console_scripts': [
            'locidex=locidex.main:main',
        ],
    },
)