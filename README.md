[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/locidex/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/locidex)
[![Conda](https://img.shields.io/conda/dn/bioconda/locidex?color=green)](https://anaconda.org/bioconda/locidex)
[![License: Apache-2.0](https://img.shields.io/github/license/phac-nml/locidex)](https://www.apache.org/licenses/LICENSE-2.0)


## Locidex
![alt text](https://github.com/phac-nml/locidex/blob/main/logo.png?raw=true)

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Quick Start](#quick-start)
- [FAQ](#faq)
- [Citation](#citation)
- [Legal](#legal)
- [Contact](#contact)

## Introduction

A common function for many tools in bacterial typing is performing similarity searching using NCBI blast. Blast provides a robust command line interface for 
constructing and using databases for similarity searching and is ubiquitous.  There are many typing applications where custom code is written around the blast
command line interface to perform searches for a variety of downstream applications. One major application within public health is identification of specific target
sequences within an assembly to perform gene-by-gene phylogenetic analysis (MLST, cgMLST, wgMLST), antimicrobial resistance gene detection, virulence gene detection, and
in silico predictions of phenotypes such as serotype. The typical approach is to bundle the search based logic with additional specialized logic for performing the desired
analysis. 

Decentralized allele calling has become a pressing concern by public health laboratories due to the increased use of whole genome sequencing (WGS) as part of outbreak detection and surveillance of 
a variety of pathogens. Gene-by-gene approaches have a variety of benefits since the introduction of 7-gene MLST approaches for species typing which include a standardized set of loci for estimating
genetic similarity between samples. This standardization allows for interoperability between different groups and also has the benefits of compression by similifying genetic comparions to use a simple
haming distance based on allele identifiers instead of a whole sequence. However, a limitation of this approach is the requirement of a centralized authority to issue unique allele identifier and this poses multiple
problems for operationalization such as privacy and connectivity. Despite this limitation PulseNet International has adopted gene-by-gene analysis as its prefered analytical approach for estimating genetic
similarity between samples for routine operations with the limitation that comparing between juristictions requires the sharing of the primary sequence data rather than the allele identifiers. 

In recent years, the concept of using cryptographic hashes of the allele sequence itself have gained traction in a variety of different allele calling software such as Chewbbaca to provide decentralized allele identifiers. Hashing the sequence yields a determinist and fixed-size 
hash value which can be compared in the same manner as integers. There are numerous hash functions with different strengths and weaknesses but MD5 digests have broad adoption in the software community and are routinely used to provide some assurance that a transferred file has arrived intact. 
The choice of md5 hash provides 16^32, possible hashes. There is a theoretical chance of hash collisions, i.e., different sequences resulting in the same hash, but as the number of viable sequences for each gene in databases is typically only tens to hundreds, this is very unlikely.

**Search**
Input Data Formats: GenBank, Fasta (of individual loci sequences)
- DNA and protein blast searches
- Md5 hashing of alleles
- Storage of results for post processing in json format

**Report**
Produce loci hash profiles in multiple formats (json, tsv, parquet)
- Filter results based on user criteria
- Multi-copy loci handling

Optional: (Not required for MVP)
Produce concatenates fasta sequences based on allele profiles

**Merge**
- Accepts list of files on command line or file of files and reads and concatenates the files into an allele profile in TSV format
- reads gz and uncompressed inputs

**Format**
- Takes common formats of gene-by-gene databases and formats them for use with locidex build module


**Build**
Builds locidex db folder structure
- Creates database configuration file
- Creates loci metadata file
- Construct blast databases (nucleotide and/or protein)


## Installation

Install the latest released version from conda:

        conda create -c bioconda -c conda-forge -n locidex locidex

Install using pip:

        pip install locidex

Install the latest master branch version directly from Github:

        pip install git+https://github.com/phac-nml/locidex.git



## Usage
If you run ``locidex``, you should see the following usage statement:



Supported input formats
=====
**GenBank**

**Fasta**

Quick start
=====


## Benchmarks

Coming soon

## FAQ

Coming soon

## Citation

Robertson, James, Wells, Matthew, Christy-Lynn, Peterson, Kyrylo Bessonov, Reimer, Aleisha, Schonfeld, Justin. LOCIDEX: Distributed allele calling engine. 2024. https://github.com/phac-nml/locidex

## Legal

Copyright Government of Canada 2023

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.


## Contact

**James Robertson**: james.robertson@phac-aspc.gc.ca
