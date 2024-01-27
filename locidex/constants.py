NT_SUB = str.maketrans('acgtrymkswhbvdnxACGTRYMKSWHBVDNX',
                       'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')

PROTEIN_ALPHA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X','*','B', 'Z','J' ]
DNA_ALPHA = ['A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']

START_CODONS = ['ATG','GTG','CTG','TTG','ATA']
STOP_CODONS = ['TAA','TAG','TTA','TCA','TGA','AGA','AGG']


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
sseq
'''.strip().split('\n')


FILE_TYPES = {
    'genbank':["gbk","genbank","gbf","gbk.gz","genbank.gz","gbf.gz"],
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
    "db_num_loci",
    "db_filenames",
]

SEARCH_RUN_DATA = {

}

DB_EXPECTED_FILES = {
    'config':'config.json',
    'meta':'meta.json',
    'nucleotide':'nucleotide.fasta',
    'protein':'protein.fasta',
}
