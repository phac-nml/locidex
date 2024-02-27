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
Blast provides a robust command line interface for  constructing and using databases for similarity searching and is ubiquitous. There are many typing 
applications where custom code is written around the blast command line interface to perform 
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
allele sequence itself have gained traction in a variety of different allele calling software such as 
[Chewbbaca](https://github.com/B-UMMI/chewBBACA) to 
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
within the database. 

Locidex is meant to be optimized for routine operation level searching where it is useful to have 
default parameters that are set for the user to have reproducibility, which is combined with flexibility to  apply 
multiple filtering parameters on the sequence store after the fact. This allows exploring different thresholds 
without the need to recompute blast searches each time and allows for different use cases of data from a common data store. 
Frequently there is a desire to include additional information about given locus such as different identifiers, 
functional properties, and phenotypic effects.  The database format of locidex allows inclusion of any number of fields 
that allow the user to describe their data which is bundled into the search result object for convenience to downstream analyses.
[Chewbbaca](https://github.com/B-UMMI/chewBBACA) is an excellent choice for an open source allele caller and provides many advanced features 
for developing, curating and using gene-by-gene schemes.  It provides a great deat of additional information regarding partial gene sequences. 
For R&D applications, this functionality can be extremely useful. However, for some operational contexts, the design of [Chewbbaca](https://github.com/B-UMMI/chewBBACA)
provides information that is not desireable and at present it has issues with multiple instances using the same database at once with 
novel allele detection enabled (https://github.com/B-UMMI/chewBBACA/issues/168). Locidex does not have the full features for a gene-by-gene software package like [Chewbbaca](https://github.com/B-UMMI/chewBBACA)
but can be used to acheive similar results while being a more generic tool kit for blast searches such as [abricate](https://github.com/tseemann/abricate)


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
    report  Filter a sequence store and produce and extract of blast results and gene profile (allele calling)
    merge   Merge a set of gene profiles into a standard "allele" profile format 
    format  Format fasta files from other MLST databases for use with locidex build
    build   Builds locidex db 


Workflows
=====
![alt text](https://github.com/phac-nml/locidex/blob/main/LocidexWorkflows.png?raw=true)

Locidex is designed to be very modular so that developers and users can mix and match different components for their individual goals.
Each tool is designed so that it can be imported as a python library to extend and implement custom behaviour. A description of each tool and 
its inputs/outputs is provided below.

**Search**

The search module is meant to use locidex formated datbase directories
Input Data Formats: GenBank, Fasta (of individual loci sequences)
- DNA and protein blast searches
- Md5 hashing of alleles 
- Storage of results for post-processing in json format

Gene annotation is notoriously inconsistent between different software and so it is **STRONGLY** recommended to use the same method
for annotation of your database and what you will use to search. ie. if using prodigal for searching, use prodigal for constructing the database.

 

    locidex search -q ./example/search/NC_003198.1.fasta -d .example/build_db_mlst_out -o ./example/search/NC_003198_fasta -n 8 --annotate
    
-run in annotation mode with a fasta input


    locidex search -q ./example/search/NC_003198.1.gbk -d .example/build_db_mlst_out -o ./example/search/NC_003198_fasta -n 8

-run with existing annotations in GenBank format

**Output**:
```
{out folder name}
├── blast
  ├── nucleotide
    ├── hsps.txt        
    └── queries.fasta
  ├── protein
    ├── hsps.txt        
    └── queries.fasta      
├── seq_store.json
└── results.json  
```
See "Sequence Store" for description of the seq_store.json output file

**Report**

Produce loci hash profiles in multiple formats (json, tsv, parquet)
- Filter results based on user criteria
- Multi-copy loci handling


    locidex report -i .example/search/seq_store.json -o ./example/report_out --name NC_003198



Optional: (Not required for MVP)
Produce concatenates fasta sequences based on allele profiles

**Merge**

Accepts list of report files on command line or file of files and reads and concatenates the files into an allele profile in TSV format (reads gz and uncompressed inputs).

        locidex merge -i ./example/merge_in/profile_1.json ./example/merge_in/profile_2.json  -o ./example/merge_out/
        
- merging multiple files provided on the command line to -i


        locidex merge -i ./example/merge_in/file_list.txt ./example/merge_out/ 

- merging a file which is a list of paths to report files

```
{out folder name}
├── profile.tsv   
└── results.json  
```

**Format**

Takes common formats of gene-by-gene databases and formats them for use with locidex build module. It accepts a directory of 
fasta files: ["fasta","fas","fa","ffn","fna","fasta.gz","fas.gz","fa.gz","ffn.gz","fna.gz"] which have the locus name of 
as the file name and allele id's are present in the fasta header separated by an underscore. ie. aroC would have the
file name aroC.fas and the header line would be >aroC_1. Additionally, format will accept a concatonated file of all 
loci in a  single fasta file which has the fasta def line as >{locus name}_{allele id}. These two formats are common with
most of the major MLST databases.

        locidex format -i ./example/format_db_mlst_in/ -o ./example/mlst_out/ 

**Output**:
```
{out folder name}
├── results.json                    
└── locidex.txt
```



**Build**

Builds locidex db folder structure
- Creates database configuration file
- Creates loci metadata file
- Construct blast databases (nucleotide and/or protein)

Takes the output of locidex format which may or may not have additional columns added. There are specific fields being looked for 
in the file which either or both are required depending on the type of db being built "dna_seq", "aa_seq". It extracts the sequence
data (nucleotide|protein) and initializes the config.json, meta.json and blast db structure which locidex search requires.

        locidex build -i ./example/build_db_mlst_in/senterica.mlst.txt -o ./example/mlst_out_db/ 

**Output**:

See - [Database structure](#Database)


## Database

Similar to  [abricate](https://github.com/tseemann/abricate), Locidex uses a fixed database structure layout. Locidex supports
nucleotide (blastn) and protein (blastp) blast searches. Locidex utilizes a few controlled fields in config.json and meta.json
but completely supports the additon of any number of additional fields that may be desired by the database builder. 
There is a nested folder structure for the blast databases which must conform to the layout below.
```
{DB folder name}
├── config.json                     #required
├── meta.json                       #required
└── blast                           #required
    └── nucleotide                  #optional but >= 1must be present
        ├──nucleotide.fasta
        ├──nucleotide.ndb
        ├──nucleotide.nhr
        ├──nucleotide.nin
        ├──nucleotide.njs
        ├──nucleotide.nsq
        ├──nucleotide.ntf
        └──nucleotide.nto 
    └──protein                      #optional but >= 1must be present
        ├── protein.fasta
        ├── protein.pdb
        ├── protein.phr
        ├── protein.pin
        ├── protein.pjs
        ├── protein.pot
        ├── protein.psq
        ├── protein.ptf
        └──protein.pto
```

**config.json**

This [JSON](https://www.json.org/json-en.html) file is responsible for encoding the metadata regarding the locidex database which will be bundled into the seq_store output
of the search module. It also determines whethere a given database is used for nucleotide and/or protein blast searches.
This data assists with provenance information regarding search results.

            {
                "db_name": "Salmonella Chewbbaca-Online cgMLST",
                "db_version": "1.0.0",
                "db_date": "2024-02-01",
                "db_author": "James Robertson",
                "db_desc": "Data obtained from: https://chewbbaca.online/species/8",
                "db_num_seqs":8558,
                "is_nucl": "True",
                "is_prot": "True",
                "nucleotide_db_name": "nucelotide",
                "protein_db_name": "protein"
            }


**meta.json**

Frequently, there is a desire to encode contextual metadata with sequence data for different analytical operations. These can
range from encoding gene and allele identifiers to phenotypic effects and beyond. The options for developers has been to encode this into the fasta
header directy with different delimeters between fields or to have an additional file with the desired contextual information separate from the sequence data.
Locidex utilizes the later approach with all of the data encoded in [JSON](https://www.json.org/json-en.html) for easy parsing by downstream applications.
This file is used by the search module to bundle this contextual information regarding the database sequences for later use without the need to pass the original
meta.json information between processes. An example stub is shown below for metadata. and any number of additional fields can be added to records.

    {
        "info": {
            "num_seqs": 8558,
            "is_cds": "True",
            "trans_table": 11,
            "dna_min_len": 72,
            "dna_max_len": 10644,
            "dna_min_ident": 80,
            "aa_min_len": 24,
            "aa_max_len": 3548,
            "aa_min_ident": 80
        },
        "meta": {
            "0": {
                "locus_name": "SAL_CGMLST_000000001",
                "locus_name_alt": "INNUENDO_cgMLST-00038301_1",
                "locus_product": NaN,
                "locus_description": NaN,
                "locus_uid": NaN,
                "dna_seq_len": 1521,
                "dna_seq_hash": "c172404ec34947c7d85f0d3e7fa3b808",
                "aa_seq_len": 507,
                "aa_seq_hash": "8a28aca12bcf1608564244a688e79db6",
                "dna_min_len": 1217,
                "dna_max_len": 1825,
                "aa_min_len": 406,
                "aa_max_len": 608,
                "dna_min_ident": 80,
                "aa_min_ident": 80
            }
        }
    }

**blast databases**

The blast folder must be present for a locidex database to be valid with nucleotide and/or protein subfolders present. 
For each folder present it must contain the fasta file used to generate the database (nucleotide.fasta|protein.fasta) and
a corresponding blast database with the folder name without ".fasta" present. The structure is designed to allow for the 
inclusion of other sequence similarity search tools at a later date with minimal modifications. Protein searching via HMMs
and other tools such as [diamond](https://github.com/bbuchfink/diamond) may be included in later releases.



Supported input formats
=====
**GenBank**

[GenBank format](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) is a widely used format for storing and sharing 
biological sequence data, primarily DNA sequences. It was developed by the National Center for Biotechnology Information 
(NCBI). GenBank format files typically have a .gb or .gbk extension. When using GenBank format as input to locidex, the
translations (if applicable) are provided by the file rather than being translated by locidex. This format is used as input
to locidex search.

The GenBank format consists of multiple componets but the relevant ones for locidex are:

Sequence Data: The complete nucleotide sequence which has been annotated.

Features: Annotations describing different regions or features of the sequence, such as coding sequences (CDS), 
regulatory elements, and other functional elements. Each feature is defined by a location (start and end positions) 
on the sequence and additional qualifiers providing details about the feature.

**Fasta**

[FASTA](https://en.wikipedia.org/wiki/FASTA_format) format is a simple and widely used text-based format for representing biological sequences, 
such as DNA, RNA, or protein sequences. It is named after the FASTA software package, which first introduced this format. 
FASTA files typically have a .fasta or .fa extension. Locidex supports both nucleotide and protein sequences as input but 
it must consist of extracted CDS, ORFs, Genes which are to be hashed. You can provide the extracted gene sequences from a tool such as [prodigal](https://github.com/hyattpd/Prodigal), [prokka](https://github.com/tseemann/prokka), 
[bakta](https://github.com/oschwengers/bakta) and an appropriate translation table if the protein coding sequences. 
You can supply contigs but unless you specify --annotate (gene annotation using pyrodigal) no extraction will occur and
the entire sequence will be treated as a single query.  This format is can be used as input to locidex search.




**Sequence store**

This is a locidex defined filetype which uses [JSON](https://www.json.org/json-en.html) to store information about the 
query sequences including their hash values along with the metadata regarding the database which was queried and its 
sequences (seq_store.json). This format is used as input to locidex report.

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

## Quick Start

MLST Example: The 7-gene MLST scheme targets from https://pubmlst.org/organisms/salmonella-spp were used as targets to extract
the full length CDS annotations from NC_003198.1 (Salmonella Typhi CT18). Note that these are not just the MLST target sequences
but the full orf and so this will differ from normal MLST results. If you want to use the subsections of the loci, you will need
to extract these using another method. These were separated into individual fasta files
for each gene but a concatonated version would also work as long as the fasta header began with the locus identifier.
Format is used to take the fasta files to create a TSV file with each of the targets and individual match thresholds for each query.
These can be modified by the user before building the database. The build function converts that TSV into a form that locidex search can use.
In order to call the alleles from the blast results of the search module, the report module is called and the samples are named based on the 
user input. If no name is specified then the base name of the file is used. Finally, the individual profiles are merged into a TSV file which 
is compatible with downstream gene profile input.

    locidex format -i ~/example/format_db_mlst_in/ -o ~/example/format_db_mlst_out/ --force
    locidex build -i ~/example/format_db_mlst_out/locidex.txt -o ~/example/build_db_mlst_out/ --force
    
    locidex search -q ~/example/search/NC_003198.1.gbk -d ~/example/build_db_mlst_out/ -o ./mlst_ncbi_annotated --force 
    locidex search -q ~/example/search/NC_003198.1.fasta -d ~/example/build_db_mlst_out/ -o ./mlst_prodigal --force --annotate
    
    locidex report -i ./mlst_ncbi_annotated/seq_store.json -o ./mlst_ncbi_annotated/report --name ncbi --force
    locidex report -i ./mlst_prodigal/seq_store.json -o ./mlst_prodigal/report --name prodigal --force
    
    locidex merge -i ./mlst_ncbi_annotated/report/profile.json .//mlst_prodigal/report/profile.json -o ./merged --force
    
    


## Benchmarks

Coming soon

## FAQ

**Can I use non-coding sequences as input to locidex?**

Yes, you can toggle off protein coding within the database config to avoid translation and protein blast searches.
However, if you want to use the "allele" identity of your locus such as an rRNA gene, then you will need to extract out the sequences
you want to match.

**Do I need to have a representitive of every allele I want to match in my database?**

No, the benefit of having dual searching with protein and dna is that you can have a sparsely populated database of sequences that 
meet your applications specific need.  This means you can deduplicate using a tool such as [cd-hit](https://sites.google.com/view/cd-hit) to remove highly similar sequences
from your database to reduce runtime and complexity. You will need to empiracly determine what level of diversity is sufficient
for your specific application.

**How does Locidex handel multicopy loci?**

Ideally a gene-by-gene scheme consists of only single copy genes but bacterial genomes are dynamic and 
genuine dulplications can occur, in addition to assembly artifacts and contamination. There are a variety of approaches
available to manages these cases. Within the 7-gene [mlst](https://github.com/tseemann/mlst) tool multiple alleles for a given locus are reported 
with a comma delimiting each allele. However, this poses an issue for calculating genetic distances since it is unclear how to treat
the multiple alleles. 1) treat the combination as a novel allele 2) blank the column 3) select the earliest allele in the database 4)
Use a similarity score to rate which is the best allele to include. The most conservative approach is to not interpret that column
by blanking it in distance calculations which results in blunting resolution which is implemented withing locidex as the conservative mode.
Alternatively, by using an approach to select only one of the loci to match will have mixed effect (options 3, 4) that can result in
inconsistencies where some isolates appeare more similar or dissimilar than they are.  The prefered method that locidex has implented as its default (normal)
mode is to combine the result into a new "allele" hash that is derived from calculating the md5 hash of the concatonated allele md5 hashes which
have been sorted alphabetically. This has the benefit of the same combination of alleles will always result in the same hash code and will match when this occurs.
Conversely, it will count a difference even when the component alleles may match between two samples.


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
