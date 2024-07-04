import pytest
import os
import locidex.classes.aligner
from locidex.constants import BlastColumns
from dataclasses import dataclass

PACKAGE_ROOT = os.path.dirname(locidex.__file__)
@dataclass
class Datasets:
    def __init__(self, tmpdir):
        self.align_inputs = {'thrA':{'thrA':{
                                'locus_name': 'thrA', 
                                'ref_id': '43', 
                                'ref_seq': 'gtgctgggccgtaatggttccgactattccgccgccgtgctggccgcctgtttacgcgctgactgctgtgaaatctggactgacgtcgatggcgtgtatacctgtgacccgcgccaggtgccggacgccagactgctgaaatcgatgtcctaccaggaagcgatggaactctcttacttcggcgccaaagtccttcaccctcgcaccataacgcctatcgcccagttccagatcccctgtctgattaaaaataccggtaatccgcaggcgccaggaacgctgatcggcgcgtccagcgacgatgataatctgccggttaaagggatctctaaccttaacaacatggcgatgtttagcgtctccggcccgggaatgaaagggatgattgggatggcggcgcgtgttttcgccgccatgtctcgcgccgggatctcggtggtgctcattacccagtcctcctctgagtacagcatcagcttctgtgtgccgcagagtgactgc', 
                                'ext_seq': 'gtgctgggccgtaatggttccgactattccgccgccgtgctggccgcctgtttacgcgctgactgctgtgaaatctggactgacgtcgatggcgtgtatacctgtgacccgcgccaggtgccggacgccaggctgttgaaatcgatgtcctaccaggaagcgatggagctctcttacttcggcgctaaagtccttcaccctcgcaccataacgcctatcgcccagttccagatcccctgtctgattaaaaataccggcaatccgcaggcgccaggaacgctgatcggcgcgtccagcgacgatgataatctgccggttaaagggatctctaaccttaacaacatggcgatgtttagcgtctccggcccgggaatgaaagggatgattgggatggcggcgcgtgttttcgccgccatgtctcgcgccgggatctcggtggtgctcattacccagtcctcctctgagtacagcatcagcttctgtgtgccgcagagtgactgc', 
                                'file': tmpdir
                            }},
                            'hemD':{'hemD': {
                                'locus_name': 'hemD', 
                                'ref_id': '16', 
                                'ref_seq': 'gcgacgttgacggaaaacgatctggtttttgccctttcacagcacgccgtcgcctttgcccacgcccaactccagcgagatggtcgaaactggcctgcgtcgccgcgctatttcgcgattggtcgcaccacggcgctcgcccttcataccgttagcgggttcgatattcgttatccattggatcgggaaatcagcgaagtcttgctacaattacctgaattacaaaatattgcgggcaaacgcgcgctgattttgcgtggcaatggcggtcgcggtcgcgaactgctgggcgaaaccctgacagctcgcggagccgaagtcagtttttgtgaatgttatcaacgaagtgcgaaacattacgatggcgcagaagaggcgatgcgctggcacactcgcggcgtaacgacgcttgttgtcaccagcggcgagatgttgcaa', 
                                'ext_seq': 'gcaacgctgacggaaaacgatctggtttttgccctttcacagcactccgtcgcctttgctcacgcccagctccagcgggatggacgaaactggcctgcgtcgccgcgctatttctcgattggccgcaccacggcgctcgcccttcataccgttagcgggttcgatattcgttatccattggatcgggaaatcagcgaagccttgctacaattacctgaattacaaaatattgcgggcaaacgcgcgctgattttgcgtggcaatggcggccgcgaactgctgggcgaaaccctgacagttcgcggagccgaagtcagtttttgtgaatgttatcaacgatgtgcgaaacattacgatggcgcggaagaagcgatgcgctggcatactcgcggcgtaacaacgcttgttgttaccagcggcgagatgttgcaa', 
                                'file': tmpdir
                            }}
                            }

@pytest.fixture
def aligner_class_fixture():
    return locidex.classes.aligner.aligner(trim_fwd=True,trim_rev=True,ext_fwd=False, ext_rev=False,fill=False, snps_only=False)

def test_perform_alignment_call_variants(tmpdir, aligner_class_fixture):
    seq_id = 'thrA'; ref_id='43'
    ext_seq_data=Datasets(tmpdir).align_inputs[seq_id]     
    result = locidex.classes.aligner.perform_alignment(ext_seq_data, tmpdir, 1)
    alignment = ext_seq_data[seq_id]['alignment']
    variants = aligner_class_fixture.call_seq(seq_id,alignment[ref_id],alignment[seq_id])
    assert list(variants.keys())==['id', 'matched_bases', 'num_diff', 'coverage', 'seq']
    assert variants['num_diff'] == 5, 'Expected value 5 for sequence differences'
    assert variants['coverage'] == 100, 'Expected value 100 for reference coverage'

    aligned_sequences_str=result[0][0] #identical to alignment variable, but needed to test parser
    aligned_sequences_dict = locidex.classes.aligner.parse_align(aligned_sequences_str)
    assert list(aligned_sequences_dict.keys()) == [ref_id, seq_id]
    assert all([True if 'alignment' in ext_seq_data[key] else False for key in ext_seq_data ]), 'Missing alignment results in locidex.classes.aligner object'
    assert len(ext_seq_data[seq_id]['alignment'].keys()) == 2, 'Should have two sequences (ref and query)'
    assert len(ext_seq_data[seq_id]['alignment'][ref_id]) == 501, 'Length should be 501 bp'
    assert len(ext_seq_data[seq_id]['alignment'][seq_id]) == 501, 'Length should be 501 bp'
    assert list(ext_seq_data[seq_id].keys()) == ['locus_name', 'ref_id', 'ref_seq', 'ext_seq', 'file', 'alignment'], 'Keys do not match'

def test_perform_seq_transform_trim_extend_gapfill(aligner_class_fixture):
    tseq = aligner_class_fixture.transform_seq( 'GTGGCAATGGCGGT',
                                                'GCGTGGCAATG---')
    assert tseq == 'GCGTGGCAATGGGT', f'Unexpected result {tseq}'
    tseq = aligner_class_fixture.trim_ends('ATGT--','ATGTAA', trim_fwd=False,trim_rev=False)
    assert 'ATGTAA' == tseq, f'Unexpected result {tseq}'
    tseq = aligner_class_fixture.trim_ends('--GTAA','ATGTAA', trim_fwd=False,trim_rev=True)
    assert 'ATGTAA' == tseq, f'Unexpected result {tseq}'
    tseq = aligner_class_fixture.trim_ends('--GT--','ATGTAA', trim_fwd=True,trim_rev=True)
    assert '--GT--' == tseq, f'Unexpected result {tseq}'
    tseq = aligner_class_fixture.extend('GCGTAA','-AGTCC',ext_fwd=False, ext_rev=False)
    assert '-AGTCC' == tseq, f'Unexpected result {tseq}'
    tseq = aligner_class_fixture.extend('GCGTAA','GAGTC-',ext_fwd=False, ext_rev=True)
    assert 'GAGTCA' == tseq, f'Unexpected result {tseq}'
    tseq = aligner_class_fixture.extend('GCGTAA','-AGTCC',ext_fwd=True, ext_rev=True)
    assert 'GAGTCC' == tseq, f'Unexpected result {tseq}'
    tseq = aligner_class_fixture.gap_fill('ATGTAA','AT--AA')
    assert 'ATGTAA' == tseq, f'Unexpected result {tseq}'

def test_coverage_and_identity(aligner_class_fixture):
    coverage = aligner_class_fixture.calc_coverage('GCGTAACCAT','GCGTAACC--')
    assert coverage == 80, "Expected different coverage value between two sequences"
    coverage = aligner_class_fixture.calc_coverage('GCGTAACCAT','-CGTA-CCA-')
    assert coverage == 70, "Expected different coverage value between two sequences"
    identity = aligner_class_fixture.count_identity('GCGTAACCAT','GCGTAACC--')
    assert identity == [8,2], "Expected different match and difference value in list"
    identity = aligner_class_fixture.count_identity('GCGTAACCAT','-CGTAACCGG')
    assert identity == [7,3], "Expected different match and difference value in list"