
from dataclasses import dataclass, asdict, fields
import pathlib
from typing import Any, Union, NamedTuple, Optional

DNA_AMBIG_CHARS = ['b', 'd', 'e', 'f', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 'u', 'v', 'w', 'x',
                   'y', 'z', '-']
DNA_IUPAC_CHARS = ['b', 'd', 'e', 'f', 'h', 'i', 'j', 'k', 'l', 'm', 'o', 'p', 'q', 'r', 's', 'u', 'v', 'w', 'x', 'y',
                   'z']

NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')

PROTEIN_ALPHA = ['a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y', 'x','*','b', 'z','j' ]
DNA_ALPHA = ['a', 'c', 'g', 't', 'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n']

START_CODONS = ['atg','gtg','ctg','ttg','ata']
STOP_CODONS = ['taa','tag','tta','tca','tga','aga','agg']


@dataclass(frozen=True)
class CharacterConstants:
    stop_codon: str = "*"


#BLAST_TABLE_COLS = '''
#qseqid
#sseqid
#qlen
#slen
#qstart
#qend
#sstart
#send
#length
#mismatch
#pident
#qcovhsp
#qcovs
#sstrand
#evalue
#bitscore
#'''.strip().split('\n')

class BlastColumns(NamedTuple):
    qseqid: str
    sseqid: str
    qlen: int
    slen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    length: int
    mismatch: str
    pident: float
    qcovhsp: float
    qcovs: float
    sstrand: str
    evalue: float
    bitscore: float

@dataclass(frozen=True)
class BlastCommands:
    # upgrading this to a string enum would be nice
    tblastn: str = "tblastn"
    blastn: str = "blastn"
    blastp: str = "blastp"

    @classmethod
    def _keys(cls) -> list:
        return [i.name for i in fields(cls)]

FILE_TYPES = {
    'genbank': [".gbk",".genbank",".gbf",".gbk.gz",".genbank.gz",".gbf.gz",".gbff",".gbff.gz"],
    'fasta': [".fasta",".fas",".fa",".ffn",".fna",".fasta.gz",".fas.gz",".fa.gz",".ffn.gz",".fna.gz"],
}


SCHEME_HEADER = [
    'uid',
    'locus_id',
    'locus_meta',
    'is_cds',
    'max_internal_stop',
    'min_dna_len',
    'max_dna_len',
    'dna_seq',
    'dna_seq_len',
    'dna_seq_hash',
    'aa_seq',
    'aa_seq_len',
    'aa_seq_hash',
    'start_codon',
    'stop_codon'
]

EXTRACT_MODES = ['snps','trim','raw','extend']

# Manifest opts for parsing


OPTION_GROUPS = {
    "db_group": ["db_name", "db_version"],
}

@dataclass
class DBConfig:
    db_name: Optional[str] = None
    db_version: Optional[str] = None
    db_date: Optional[str] = None
    db_author: Optional[str] = None
    db_desc: Optional[str] = None
    db_num_seqs: Optional[Union[str, int]] = None
    is_nucl: Optional[bool] = None
    is_prot: Optional[bool] = None
    nucleotide_db_name: Optional[str] = None
    protein_db_name: Optional[str] = None

    def __getitem__(self, name: str) -> Any:
        return getattr(self, str(name))
    
    def __setitem__(self, key: str, value: str) -> None:
        setattr(self, key, value)
    
    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def _keys(cls) -> list:
        return [i.name for i in fields(cls)]

    def keys(self) -> list:
        return [i.name for i in fields(self)]


@dataclass(frozen=True)
class DBFiles:
    meta_file: str = "meta.json"
    config_file: str = "config.json"
    results_file: str = "results.json"
    blast_dir: str = "blast"

    @classmethod
    def items(cls):
        return [(i.name, pathlib.Path(i.default)) for i in fields(cls)]

@dataclass(frozen=True)
class ManifestFields:
    db_path: str = "path"
    config_data: str = "config"

SEARCH_RUN_DATA = {

}


DB_EXPECTED_FILES = {
    'config':'config.json',
    'meta':'meta.json',
    'nucleotide':'nucleotide.fasta',
    'protein':'protein.fasta',
}

@dataclass
class MetadataFields:
    num_seqs: int
    is_cds: bool
    trans_table: int
    dna_min_len: int
    dna_max_len: int
    dna_min_ident: float
    aa_min_len: int
    aa_max_len: int
    aa_min_ident: float
    
    def to_dict(self):
        return asdict(self)



class LocidexDBHeader(NamedTuple):
    seq_id: str
    locus_name: str
    locus_name_alt: str 
    locus_product: str
    locus_description: str
    locus_uid: str
    dna_seq: str
    dna_seq_len: int
    dna_seq_hash: str
    aa_seq: Optional[str]
    aa_seq_len: Optional[int]
    aa_seq_hash: Optional[str]
    dna_min_len: int
    dna_max_len: int
    aa_min_len: Optional[int]
    aa_max_len: Optional[int]
    dna_min_ident: float
    aa_min_ident: Optional[float]
    min_dna_match_cov: int
    min_aa_match_cov: Optional[int]
    count_int_stops: int
    dna_ambig_count: int

