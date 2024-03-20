[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/locidex/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/locidex)
[![Conda](https://img.shields.io/conda/dn/bioconda/locidex?color=green)](https://anaconda.org/bioconda/locidex)
[![License: Apache-2.0](https://img.shields.io/github/license/phac-nml/locidex)](https://www.apache.org/licenses/LICENSE-2.0)

<img src="https://github.com/phac-nml/locidex/blob/dev/assets/logo.png" width = "180" height="140">

# Locidex
# Introduction
A common function for many tools in bacterial typing is performing similarity searching using NCBI [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi). Blast provides a robust command line interface for constructing and using databases for similarity searching and is ubiquitous. There are many typing applications where custom code is written around the blast command line interface to perform searches for a variety of downstream applications. For instance, identification of specific target sequences within an assembly to perform gene-by-gene phylogenetic analysis (MLST, cgMLST, wgMLST), antimicrobial resistance gene detection, virulence gene detection, and in silico predictions of phenotypes such as serotype is a major application within public health. The typical approach is to bundle the search-based logic with additional specialized logic for performing the desired analysis.

Decentralized allele calling has become a pressing concern by public health laboratories due to the increased use of whole genome sequencing (WGS) as part of outbreak detection and surveillance of a variety of pathogens. Gene-by-gene approaches have a variety of benefits for species typing which include a standardized set of loci for estimating genetic similarity between samples. This standardization allows for interoperability between different groups and also has the benefit of compression, simplifying genetic comparisons to use a simple hamming distance based on allele identifiers instead of a whole sequence. However, a limitation of this approach is the requirement of a centralized authority to issue unique allele identifiers and this poses multiple problems for operationalization such as privacy and connectivity. Despite this limitation PulseNet International has adopted gene-by-gene analysis as its preferred analytical approach for estimating genetic similarity between samples for routine operations with the limitation that comparing between jurisdictions requires the sharing of the primary sequence data rather than the allele identifiers.

In recent years, the concept of using cryptographic hashes of the allele sequence itself have gained traction in a variety of different allele calling software, such as [Chewbbaca](https://github.com/B-UMMI/chewBBACA), to provide decentralized allele identifiers. Hashing the sequence yields a determinist and fixed-size hash value which can be compared in the same manner as integers. There are numerous hash functions with different strengths and weaknesses but MD5 digests have broad adoption in the software community and are routinely used to provide some assurance that a transferred file has arrived intact. The choice of md5 hash provides 16^32, possible hashes. There is a theoretical chance of hash collisions, i.e., different sequences resulting in the same hash, but as the number of allele sequences for each gene in databases is relatively low, this should be an uncommon occurrence. Collisions in this case would result in profiles appearing more similar than they truly are at the sequence level. In addition, the chances of multiple occurrences of collisions within a profile would be infinitely small.

The motivation for developing locidex is the need a common searching engine for various loci based typing applications such as: gene-by-gene (mlst, cgMLST, wgMLST, rmlst), in silico serotyping, gene-based phenotype predictions (amr, virulence, pathotype, toxin typing), marker-based typing (16S). The tool must provide custom criteria filtering by loci, and produce multiple formats for downstream applications. It must be compatible with an HSP environment and not encounter any locking issues where multiple processes may try to change the data at the same time. [THIS SECTION WILL NEED EDITING]The logic for allele calling is greatly simplified by leveraging existing annotations from tools such as [prodigal](https://github.com/hyattpd/Prodigal), [prokka](https://github.com/tseemann/prokka), [bakta](https://github.com/oschwengers/bakta) to delineate the boundaries of the sequences to be queried and hashed to produce allele identifiers. A common issue in matching applications is that ranges of identity and coverage for a match will vary by locus and so locidex builds into its database structure control over these attributes at a locus level allowing for high variability databases to be used without building custom logic downstream. This is particularly important when lengths of loci can exhibit considerable variability as is the case for genes of interest for typing applications. This provides greater flexibility for the designation of ideal thresholds for a given application. However, these values can be overridden using the report module filtering parameters as well as by modifying the values within the database. [END]

[Chewbbaca](https://github.com/B-UMMI/chewBBACA) is an excellent choice for an open source allele caller and provides many advanced features for developing, curating and using gene-by-gene schemes. It provides a great deat of additional information regarding partial gene sequences. For R&D applications, this functionality can be extremely useful. However, for some operational contexts, the design of [Chewbbaca](https://github.com/B-UMMI/chewBBACA) provides undesirable information and at present it has issues with multiple instances using the same database at once with novel allele detection enabled ([B-UMMI/chewBBACA#168](https://github.com/B-UMMI/chewBBACA/issues/168)). Locidex is meant to be optimized for routine operation level searching where it is useful to have default parameters that are set for the user to have reproducibility combined with flexibility to apply multiple filtering parameters on the sequence store after the fact. This allows exploring different thresholds without the need to recompute blast searches. In addition, there is often a desire to include additional information about a given locus such as different identifiers, functional properties, and phenotypic effects. The database format of locidex allows inclusion of any number of fields bundled into a search result object for users to describe their data conveniently during downstream analysis. This functionality allows for different use cases of data from a common data store. Locidex does not have the full features for a gene-by-gene software package like [Chewbbaca](https://github.com/B-UMMI/chewBBACA) but can be used to achieve similar results while being a more generic tool kit for blast searches, similar to [abricate](https://github.com/tseemann/abricate).

## Citation

Robertson, James, Wells, Matthew, Christy-Lynn, Peterson, Kyrylo Bessonov, Reimer, Aleisha, Schonfeld, Justin. LOCIDEX: Distributed allele calling engine. 2024. [https://github.com/phac-nml/locidex](https://github.com/phac-nml/locidex)
## Contact

For any questions, issues or comments please make a Github issue or reach out to [**James Robertson**](james.robertson@phac-aspc.gc.ca).

# Install

Install the latest released version from conda:

        conda create -c bioconda -c conda-forge -n locidex locidex

Install using pip:

        pip install locidex

Install the latest master branch version directly from Github:

        pip install git+https://github.com/phac-nml/locidex.git

### Compatibility

[List out Dependencies and/or packages as appropriate]

# Getting Started

## Usage

		locidex <command> [options] <required arguments>

### Commands

Locidex uses the following commands:

1. **search** - query a set of ORFs, and genes against a database to produce a sequence store for downstream processing
2. **extract** - extract loci from a genome based on a locidex database
3. **report** - filter a sequence store and produce an extract of blast results and gene profile (allele calling)
4. **merge** - merge a set of gene profiles into a standard 'allele' profile format
5. **format** - format fasta files from other MLST databases for use with locidex build
6. **build** - builds a locidex databse

## Configuration and Settings:

Locidex is designed to be very modular so that developers and users can mix and match different components for their individual goals. Each tool is designed so that it can be imported as a python library to extend and implement custom behaviour. A description of each tool and its inputs/outputs is provided below.

The below figure shows a general workflow for each of the locidex commands:
![img1](https://github.com/phac-nml/locidex/blob/dev/assets/locidex_workflow_mermaid-20240318.png)

### Search

The search module is meant to use locidex formatted database directories. 

- DNA and protein blast searches
- Md5 hashing of alleles
- Storage of results for post-processing in json format

Gene annotation is notoriously inconsistent between different software and so it is **STRONGLY** recommended to use the same method of annotation for your database and what you will use to search. ie. if using prodigal for searching, use prodigal for constructing the database.

EXAMPLE: Run search in annotation mode with a fasta input:

		locidex search -q ./example/search/NC_003198.1.fasta -d ./example/build_db_mlst_out -o ./example/search/NC_003198_fasta -n 8 --annotate

EXAMPLE: Run search with existing annotations in GenBank format:

		locidex search -q ./example/search/NC_003198.1.gbk -d ./example/build_db_mlst_out -o ./example/search/NC_003198_fasta -n 8

#### Input

Accepted input Data Formats: GenBank, Fasta (of individual loci sequences)

#### Output:

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

### Extract

The extract module is meant to use locidex formatted database directories to get sequences of individual loci based on a locidex formatted database. The extract module operates in four different modes:

1) raw: sequences are directly extracted from the assembly with no further processing.
2) trim: any leading or trailing bases which are not present in the db match are trimmed from the sequence. 
3) snp: This will apply only nucleotide variants to the reference allele which can be very useful for nanopore assemblies where indels are common and unlikely to be real. 
4) extend : This mode will fill in any terminal sequence missing from the sequence based on the matched reference allele.

> [!Note]
> Modes involving processing (trim, snp and extend) involve pairwise mafft alignment of the extracted sequence with its best blast hit in the database.

![extract command options](https://github.com/phac-nml/locidex/blob/main/assets/LocidexExtract.png?raw=true)

#### Input

Input Data Formats: Fasta (contigs)

Gene annotation is notoriously inconsistent between different software, and so we implemented the extract module to enable consistent selection of loci sequences from an input genome.

EXAMPLE: to extract loci sequences from an input genome, reporting just extracted sequences and skipping any post processing (mode=`raw`)

		locidex extract --mode raw -i ./example/search/NC_003198.1.fasta -d .example/build_db_mlst_out -o ./example/search/NC_003198_fasta -n 8 

#### Output

```
{out folder name}
├── blast
    ├── hsps.txt        
├── blast_db
    ├── contigs.fasta.ndb  
    ├── contigs.fasta.nhr
    ├── contigs.fasta.nin 
    ├── contigs.fasta.njs
    ├── contigs.fasta.not
    ├── contigs.fasta.nsq
    ├── contigs.fasta.ntf
    └──contigs.fasta.nto
├── filtered.hsps.txt
├── processed.extracted.seqs.fasta #optional sequences with trimming, gapp filling an snp only based on options selected
├── raw.extracted.seqs.fasta #exact extracted sequences
└── results.json  
```

### Report

Produce loci hash profiles in multiple formats (json, tsv, parquet)

- Filter results based on user criteria
- Multi-copy loci handling

**Optional:** (Not required for MVP) Produce concatenated fasta sequences based on allele profiles
#### Input

A Sequence store (`seq_store.json`) object produced by the 'search' function.

	    locidex report -i .example/search/seq_store.json -o ./example/report_out --name NC_003198

#### Output

[INSERT REPORT OUTPUT]

### Merge

Reads and concatenates report files into an allele profile in TSV format.

#### Input

Can list report files on command line or provide a 'file of files' (FOF). ('gz' compressed and uncompressed files are excepted)

EXAMPLE: merging multiple files provided on the command line to -i

		locidex merge -i ./example/merge_in/profile_1.json ./example/merge_in/profile_2.json  -o ./example/merge_out/

EXAMPLE: merging files provided through a list of paths to report files

		  locidex merge -i ./example/merge_in/file_list.txt ./example/merge_out/ 

#### Output

```
{out folder name}
├── profile.tsv   
└── results.json  
```

### Format

Takes common formats of gene-by-gene databases and formats them for use with locidex build module.
#### Input

Accepts two formats common with most of the major MLST databases:

1. a directory of fasta files: ["fasta","fas","fa","ffn","fna","fasta.gz","fas.gz","fa.gz","ffn.gz","fna.gz"] with "locus name" as the file name and allele id's are present in the fasta header separated by an underscore. ie. aroC would have the file name aroC.fas and the header line would be >aroC_1. 
2. a concatonated file of all loci in a single fasta file which has the fasta def line as `>{locus name}_{allele id}`. These two formats are common with most of the major MLST databases.

	    locidex format -i ./example/format_db_mlst_in/ -o ./example/mlst_out/ 

#### Output

```
{out folder name}
├── results.json                    
└── locidex.txt
```

### Build

Builds locidex db folder structure

- Creates database configuration file
- Creates loci metadata file
- Construct blast databases (nucleotide and/or protein)

#### Input

Takes the output of **locidex format** (may or may not have additional columns added). There are specific fields being looked for in the file which either or both are required depending on the type of db being built "dna_seq", "aa_seq".
[I THINK THIS NEEDS LOOKING AT]


		locidex build -i ./example/build_db_mlst_in/senterica.mlst.txt -o ./example/mlst_out_db/ 

#### Output

This command extracts the sequence data (nucleotide|protein) and initializes the config.json, meta.json and blast db structure which locidex search requires.
See - [Database structure](/README.md#Database) for further information.

## Example workflow

MLST Example: The 7-gene MLST scheme targets from [https://pubmlst.org/organisms/salmonella-spp](https://pubmlst.org/organisms/salmonella-spp) were used as targets to extract the full length CDS annotations from NC_003198.1 (Salmonella Typhi CT18). Sequences were separated into individual fasta files for each gene, though a concatonated version would also work as long as the fasta header began with the locus identifier. 

> [!Note]
> The extracted  CDS annotations are not just the MLST target sequences but the full orf and so this will differ from normal MLST results. If you want to use the traditional subsections of each loci, you will need to extract these using another method. 

`locidex format` is used to create a TSV file containing the sequence of each of the targets and individual match thresholds for each query. These can be modified by the user before building the database. 

		locidex format -i ~/example/format_db_mlst_in/ -o ~/example/format_db_mlst_out/ --force

The `locidex build` converts that TSV into a form that `locidex search` can use.

		locidex build -i ~/example/format_db_mlst_out/locidex.txt -o ~/example/build_db_mlst_out/ --force

`locidex search` is used to query against the database to produce a sequence store (two examples are provided here to show the use of genbank annotations or prodigal results).

		locidex search -q ~/example/search/NC_003198.1.gbk -d ~/example/build_db_mlst_out/ -o ./mlst_ncbi_annotated --force 

		locidex search -q ~/example/search/NC_003198.1.fasta -d ~/example/build_db_mlst_out/ -o ./mlst_prodigal --force --annotate

To call the alleles from the blast results of the search module, the `locidex report` is called and the samples are named based on the user input. If no name is specified then the base name of the file is used. (two examples are provided to show the difference when using genbank annotations vs prodigal results)

		locidex report -i ./mlst_ncbi_annotated/seq_store.json -o ./mlst_ncbi_annotated/report --name ncbi --force

		locidex report -i ./mlst_prodigal/seq_store.json -o ./mlst_prodigal/report --name prodigal --force

Finally, the individual profiles are merged into a TSV file which is compatible with downstream gene profile input via `locidex merge`

		locidex merge -i ./mlst_ncbi_annotated/report/profile.json .//mlst_prodigal/report/profile.json -o ./merged --force

## Database structure

Similar to [abricate](https://github.com/tseemann/abricate), Locidex uses a fixed database structure layout. Locidex supports nucleotide (blastn) and protein (blastp) blast searches. It utilizes a few controlled fields in config.json and meta.json but completely supports the additon of any number of additional fields that may be desired by the database builder. There is a nested folder structure for the blast databases which must conform to the layout below.

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

### config.json

This [JSON](https://www.json.org/json-en.html) file is responsible for encoding the metadata regarding the locidex database which will be bundled into the seq_store output of the search module. It also determines whether a given database is used for nucleotide and/or protein blast searches. This data assists with provenance information regarding search results.

```
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
```

### meta.json

Frequently, there is a desire to encode contextual metadata with sequence data for different analytical operations. These can range from encoding gene and allele identifiers to phenotypic effects and beyond. The options for developers has been to encode this into the fasta header directy with different delimeters between fields or to have an additional file with the desired contextual information separate from the sequence data. Locidex utilizes the later approach with all of the data encoded in [JSON](https://www.json.org/json-en.html) for easy parsing by downstream applications. This file is used by the search module to bundle this contextual information regarding the database sequences for later use without the need to pass the original meta.json information between processes. An example stub is shown below for metadata. and any number of additional fields can be added to records.

```
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
```

### blast datbases

The blast folder must be present for a locidex database to be valid with nucleotide and/or protein subfolders present. For each folder present, it must contain the fasta file used to generate the database (nucleotide.fasta|protein.fasta) and a corresponding blast database with the folder name without ".fasta" present. The structure is designed to allow for the inclusion of other sequence similarity search tools at a later date with minimal modifications. Protein searching via HMMs and other tools such as [diamond](https://github.com/bbuchfink/diamond) may be included in later releases.

# Troubleshooting and FAQs:

## FAQ

**Can I use non-coding sequences as input to locidex?**

Yes, you can toggle off protein coding within the database config to avoid translation and protein blast searches. However, if you want to use the "allele" identity of your locus such as an rRNA gene, then you will need to extract out the sequences you want to match.

**Do I need to have a representitive of every allele I want to match in my database?**

No, the benefit of having dual searching with protein and dna is that you can have a sparsely populated database of sequences that meet your applications specific need. This means you can deduplicate using a tool such as [cd-hit](https://sites.google.com/view/cd-hit) to remove highly similar sequences from your database to reduce runtime and complexity. You will need to empiracly determine what level of diversity is sufficient for your specific application.

**How does Locidex handel multicopy loci?**

Ideally a gene-by-gene scheme consists of only single copy genes but bacterial genomes are dynamic and genuine dulplications can occur, in addition to assembly artifacts and contamination. There are a variety of approaches available to manages these cases. Within the 7-gene [mlst](https://github.com/tseemann/mlst) tool multiple alleles for a given locus are reported with a comma delimiting each allele. However, this poses an issue for calculating genetic distances since it is unclear how to treat the multiple alleles. There are several common methods for how to treat multiple alleles:

1) treat the combination as a novel allele 
2) blank the column 
3) select the earliest allele in the database 
4) Use a similarity score to rate which is the best allele to include. 

The most conservative approach is to not interpret that column by blanking it in distance calculations which results in blunting resolution which is implemented within locidex as the conservative mode. Alternatively, by using an approach to select only one of the loci to match will have mixed effects (options 3, 4) that can result in inconsistencies where some isolates appear more similar or dissimilar than they are. The preferred method that locidex has implemented as its [DEFAULT] (normal) mode is to combine the result into a new "allele" hash that is derived from calculating the md5 hash of the concatenated allele md5 hashes, sorted alphabetically. This has the benefit of the same combination of alleles  resulting in the same hash code and will match when this occurs. Conversely, it will count a difference even when individual component alleles may match between two samples.

# Benchmarking

Coming soon.

# Legal and Compliance Information:

Copyright Government of Canada 2023

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# Updates and Release Notes:
Release notes highlighting new features and changes.
