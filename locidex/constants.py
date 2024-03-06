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


BLAST_TABLE_COLS = '''
qseqid
sseqid
qlen
slen
qstart
qend
sstart
send
length
mismatch
pident
qcovhsp
qcovs
sstrand
evalue
bitscore
'''.strip().split('\n')


FILE_TYPES = {
    'genbank':["gbk","genbank","gbf","gbk.gz","genbank.gz","gbf.gz","gbff","gbff.gz"],
    'fasta':["fasta","fas","fa","ffn","fna","fasta.gz","fas.gz","fa.gz","ffn.gz","fna.gz"],
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

DB_CONFIG_FIELDS = [
    "db_name",
    "db_version",
    "db_date",
    "db_author",
    "db_desc",
    "db_num_seqs",
    "is_nucl",
    "is_prot",
    "nucleotide_db_name",
    "protein_db_name",
]

SEARCH_RUN_DATA = {

}

FORMAT_RUN_DATA = {

}

DB_EXPECTED_FILES = {
    'config':'config.json',
    'meta':'meta.json',
    'nucleotide':'nucleotide.fasta',
    'protein':'protein.fasta',
}

LOCIDEX_DB_HEADER = [
    'seq_id',
    'locus_name',
    'locus_name_alt',
    'locus_product',
    'locus_description',
    'locus_uid',
    'dna_seq',
    'dna_seq_len',
    'dna_seq_hash',
    'aa_seq',
    'aa_seq_len',
    'aa_seq_hash',
    'dna_min_len',
    'dna_max_len',
    'aa_min_len',
    'aa_max_len',
    'dna_min_ident',
    'aa_min_ident',
    'min_dna_match_cov',
    'min_aa_match_cov',
    'count_int_stops',
    'dna_ambig_count'

]
