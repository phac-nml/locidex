import sys

from locidex.classes.mafft import mafft
from multiprocessing import Pool, cpu_count
import os
from locidex.classes import run_command


def stage_data(seq_data,work_dir):
    d = 0
    for id in seq_data:
        record = seq_data[id]
        locus_name = record['locus_name']
        ref_id = record['ref_id']
        ref_seq = record['ref_seq']
        ext_seq = record['ext_seq']
        out_file = os.path.join(work_dir, f"{d}-{locus_name}.fas")
        write_seq({ref_id: ref_seq, id: ext_seq}, out_file)
        seq_data[id]['file'] = out_file
        d += 1


def write_seq(data, out_file):
    with open(out_file, 'w') as oh:
        for id in data:
            oh.write(">{}\n{}\n".format(id, data[id]))



def align(input_fasta,params={'globalpair':''}):
    command = ['mafft']
    for p in params:
        command.append(f'--{p}')
        command.append(params[p])
    command.append(input_fasta)
    command = " ".join([str(x) for x in command])
    (stdout, stderr) = run_command(command)
    return (stdout,stderr)


def parse_align(align):
    align = align.split('\n')
    seqs = {}
    for row in align:
        if len(row) == 0:
            continue
        if row[0] == '>':
            id = row[1:]
            seqs[id] = []
            continue
        seqs[id].append(row.lower())
    for id in seqs:
        seqs[id] = "".join([str(x) for x in seqs[id]])
    return seqs

def perform_alignment(seq_data,work_dir,num_threads=1):
    try:
        sys_num_cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        sys_num_cpus = cpu_count()

    if num_threads > sys_num_cpus:
        num_threads = sys_num_cpus

    stage_data(seq_data, work_dir)

    pool = Pool(processes=num_threads)

    results = []
    for id in seq_data:
        record = seq_data[id]
        results.append(pool.apply_async(align, args=((record['file'],))))

    pool.close()
    pool.join()

    r = []
    for x in results:
        if isinstance(x, dict):
            r.append(x)
        else:
            r.append(x.get())

    ids = list(seq_data.keys())
    for i in range(0, len(r)):
        id = ids[i]
        if len(r[i]) == 0:
            continue
        seq_data[id]['alignment'] = parse_align(r[i][0])

    return r


class aligner:
    def __init__(self,trim_fwd=True,trim_rev=True,ext_fwd=True, ext_rev=True,fill=True,snps_only=False):
        self.trim_fwd = trim_fwd
        self.trim_rev = trim_rev
        self.ext_fwd = ext_fwd
        self.ext_rev = ext_rev
        self.fill = fill
        self.snps_only = snps_only


    def transform_seq(self,seq1,seq2,trim_fwd=True,trim_rev=True,ext_fwd=True, ext_rev=True,gap_fill=True,snps_only=False):
        if snps_only:
            length = len(seq1)
            tseq = list(seq1)
            for i in range(0,length):
                b1 = seq1[i]
                b2 = seq2[i]
                if b1 != '-' and b2 != '-':
                    tseq[i] = b2

            tseq = "".join([str(x) for x in tseq])
        else:
            tseq = self.trim_ends(seq1,seq2,trim_fwd,trim_rev)
            tseq = self.extend(seq1, tseq, ext_fwd, ext_rev)
            if gap_fill:
                tseq = self.gap_fill(seq1, tseq)
        return tseq


    def trim_ends(self,seq1,seq2,trim_fwd=True,trim_rev=True):
        seq = list(seq2)
        is_end = False
        length = len(seq1)
        if trim_fwd:
            for i in range(0,length):
                b1 = seq1[i]
                if b1 == '-':
                    is_end = True
                if is_end and b1 == '-':
                    seq[i] = '-'
                else:
                    is_end = False
                    break

        is_end = False
        length = len(seq1)
        if trim_rev:
            for i in reversed(range(0,length)):
                b1 = seq1[i]
                if b1 == '-':
                    is_end = True
                if is_end and b1 == '-':
                    seq[i] = '-'
                else:
                    break
        return "".join([str(x) for x in seq])



    def extend(self, seq1, seq2,ext_fwd=True, ext_rev=True):
        seq = list(seq2)
        is_end = True
        length = len(seq1)

        if ext_fwd:
            for i in range(0, length):
                b1 = seq1[i]
                b2 = seq2[i]
                if b1 == '-':
                    is_end = False
                if is_end and b1 != '-' and b2 == '-':
                    seq[i] = b1


        is_end = True
        length = len(seq1)
        if ext_rev:
            for i in reversed(range(0, length)):
                b1 = seq1[i]
                b2 = seq2[i]
                if b1 == '-':
                    is_end = False

                if is_end and b1 != '-'  and b2 == '-':
                    seq[i] = b1

        return "".join([str(x) for x in seq])

    def gap_fill(self, seq1, seq2):
        seq = []
        is_end = True
        length = len(seq1)
        for i in range(0, length):
            b1 = seq1[i]
            b2 = seq2[i]
            if b1 == '-':
                is_end = False
            if is_end and b2 == '-':
                b2 = b1
            seq.append(b2)
        return "".join([str(x) for x in seq])

    def calc_coverage(self,seq1,seq2):
        coverage = 0
        total = 0
        for idx,base in enumerate(seq1):
            b1 = seq1[idx]
            b2 = seq2[idx]
            if b1 != '-':
                total+=1
                if b2 != '-':
                    coverage+=1
        if total > 0:
            return coverage / total * 100
        return 0


    def call_seq(self,id, ref_seq, query_seq):
        seq_record = {
            'id': id,
            'matched_bases': 0,
            'num_diff': 0,
            'coverage': 0,
            'seq': ''
        }
        if self.trim_fwd or self.trim_rev or self.ext_fwd or self.ext_rev or self.gap_fill:
            query_seq = self.transform_seq(ref_seq,query_seq,self.trim_fwd,self.trim_rev,
                                      self.ext_fwd, self.ext_rev,self.fill,self.snps_only)
        (match, diff ) = self.count_identity(ref_seq,query_seq)
        seq_record['matched_bases'] = match
        seq_record['num_diff'] = diff
        seq_record['seq'] = query_seq
        return seq_record

    def count_identity(self,seq1,seq2):
        match = 0
        diff = 0
        for idx, base in enumerate(seq1):
            if seq1[idx] != '-' and seq1[idx] == seq2[idx]:
                match+=1
            else:
                diff+=1
        return [match, diff]
