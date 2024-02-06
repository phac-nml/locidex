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

A common function for many tools in bacterial typing is performing similarity searching using NCBI [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi). 
Blast provides a robust command line interface for  constructing and using databases for similarity searching and is ubiquitous.  
There are many typing applications where custom code is written around the blast command line interface to perform 
searches for a variety of downstream applications. One major application within public health is identification of 
specific target sequences within an assembly to perform gene-by-gene phylogenetic analysis 
(MLST, cgMLST, wgMLST), antimicrobial resistance gene detection, virulence gene detection, and in silico predictions of 
phenotypes such as serotype. The typical approach is to bundle the search-based logic with additional specialized logic 
for performing the desired analysis.   

Decentralized allele calling has become a pressing concern by public health laboratories due to the increased use of 
whole genome sequencing (WGS) as part of outbreak detection and surveillance of  a variety of pathogens. 
Gene-by-gene approaches have a variety of benefits since the introduction of 7-gene MLST approaches for species typing 
which include a standardized set of loci for estimating genetic similarity between samples. This standardization allows 
for interoperability between different groups and also has the benefits of compression by simplifying genetic comparisons 
to use a simple hamming distance based on allele identifiers instead of a whole sequence. However, a limitation of this 
approach is the requirement of a centralized authority to issue unique allele identifier and this poses multiple problems 
for operationalization such as privacy and connectivity. Despite this limitation PulseNet International has adopted 
gene-by-gene analysis as its preferred analytical approach for estimating genetic similarity between samples for 
routine operations with the limitation that comparing between jurisdictions requires the sharing of the primary 
sequence data rather than the allele identifiers.   

In recent years, the concept of using cryptographic hashes of the 
allele sequence itself have gained traction in a variety of different allele calling software such as Chewbbaca to 
provide decentralized allele identifiers. Hashing the sequence yields a determinist and fixed-size hash value which 
can be compared in the same manner as integers. There are numerous hash functions with different strengths and weaknesses 
but MD5 digests have broad adoption in the software community and are routinely used to provide some assurance that a 
transferred file has arrived intact.  The choice of md5 hash provides 16^32, possible hashes. 
There is a theoretical chance of hash collisions, i.e., different sequences resulting in the same hash, but as the number 
of allele sequences for each gene in databases is relatively low, this should be an uncommon occurrence. Collisions in 
this case would just have the consequence of having a profile appear more similar than they truly are at the sequence 
level but the chances of multiple occurrences of collisions within a profile would be infinitely small.  

The motivation for developing locidex is the need a common searching engine for various loci based typing applications 
such as: gene-by-gene (mlst, cgMLST, wgMLST, rmlst), in silico serotyping, gene-based phenotype predictions 
(amr, virulence, pathotype, toxin typing), marker-based typing (16S).  The tool needs to provide custom criteria filtering 
by loci, and ability to produce multiple formats for different downstream applications to use. The tool needs to be 
compatible with an HSP environment and not encounter any locking issues where multiple processes may try to change the 
data at the same time.  The logic for allele calling is greatly simplified by leveraging existing annotations from tools 
such as [prodigal](https://github.com/hyattpd/Prodigal), [prokka](https://github.com/tseemann/prokka), 
[bakta](https://github.com/oschwengers/bakta) to delinate the boundaries of the sequences to be queried and hashed to produce allele 
identifiers. A common issue in matching applications is that ranges of identity and coverage for a match will vary by locus 
and so locidex builds into its database structure control over these attributes at a locus level which allows for 
high variability databases to be used without building custom logic downstream. This is particularly important when lengths 
of loci can exhibit considerable variability as is the case for genes of interest for typing applications. 
This provides greater flexibility for the designation of ideal thresholds for a given application. 
However, these values can be overridden using the report module filtering parameters as well as by modifying the values 
within the database. Locidex is meant to be optimized for routine operation level searching where it is useful to have 
default parameters that are set for the user to have reproducibility, which is combined with flexibility to  apply 
multiple filtering parameters on the sequence store after the fact. This allows exploring different thresholds 
without the need to recompute blast searches each time and allows for different use cases of data from a common data store. 
Frequently there is a desire to include additional information about given locus such as different identifiers, 
functional properties, and phenotypic effects.  The database format of locidex allows inclusion of any number of fields 
that allow the user to describe their data which is bundled into the search result object for convenience to downstream analyses.  

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

    Usage: locidex <command> [options] <required arguments>
    
    To get minimal usage for a command use:
    gas command
    
    To get full help for a command use one of:
    locidex command -h
    locidex command --help
    
    
    Available commands:
    
    search  Query set of ORFs, Genes against a database to produce a sequence store for downstream processing
    report  Filter a sequence store and produce and extract of blast results and gene profile
    merge   Merge a set of gene profiles into a standard profile format
    format  Not Implemented
    build   Not Implemented



Supported input formats
=====
**GenBank**

[GenBank format](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) is a widely used format for storing and sharing 
biological sequence data, primarily DNA sequences. It was developed by the National Center for Biotechnology Information 
(NCBI). GenBank format files typically have a .gb or .gbk extension. When using GenBank format as input to locidex, the
translations (if applicable) are provided by the file rather than being translated by locidex.

The GenBank format consists of multiple componets but the relevant ones for locidex are:

Sequence Data: The complete nucleotide sequence which has been annotated.

Features: Annotations describing different regions or features of the sequence, such as coding sequences (CDS), 
regulatory elements, and other functional elements. Each feature is defined by a location (start and end positions) 
on the sequence and additional qualifiers providing details about the feature.

**Fasta**

[FASTA](https://en.wikipedia.org/wiki/FASTA_format) format is a simple and widely used text-based format for representing biological sequences, 
such as DNA, RNA, or protein sequences. It is named after the FASTA software package, which first introduced this format. 
FASTA files typically have a .fasta or .fa extension. Locidex supports both nucleotide and protein sequences as input but 
it must consist of extracted CDS, ORFs, Genes which are to be hashed. You can supply contigs but no extraction will occur and
the entire sequence will be treated as a single query. You will need to provide the extracted gene sequences from a tool such as [prodigal](https://github.com/hyattpd/Prodigal), [prokka](https://github.com/tseemann/prokka), 
[bakta](https://github.com/oschwengers/bakta) and an appropriate translation table if the protein coding sequences.

**Sequence store**

This is a locidex defined filetype which uses [JSON](https://www.json.org/json-en.html) to store information about the 
query sequences including their hash values along with the metadata regarding the database which was queried and its 
sequences. This format is used as input to locidex report.

    {
            'db_info': { database configuration file data},
            'db_seq_info': { database sequence metadata },
            'query_data': {
                'sample_name':'User supplied or query filename',
                'query_seq_data': { hash values and extracted data on query sequences},
                'query_hit_columns': [names of blast fields used in search],
                'query_hits': { blast hit data indexed by query sequence },
                "locus_profile":{ assignment of query sequences to database loci }
            }
    }


**Gene Profile**

This is a locidex defined filetype which uses [JSON](https://www.json.org/json-en.html) to store information about the 
query sequence hashes that are associated to a locus. All loci in the database will appear in database order
and are included whether they have a query match or not. This format is used as input to locidex merge.

    sample_name : {
                'locus_1': '219699eecfd1176d6c6d5409b60ed556',
                'locus_2': '758146bd697a858c8cdca7322b542881',
                'locus_3': '',                                          #Missing locus
                'locus_4': 'e20aacb7efe5289a5626377a45ddcfa7',
                'locus_5': '9e5bbbdb0ab2f5ec85088bacd0167b85',
                'locus_6': '272e84a4ce4d8a376fd5a41bec94ad36',
                'locus_7': 'cc1311359fd1cd1803722ecffec279fc,50605327a347ca7d540722e69bf2204f', #Multiple queries match locus
    }



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
