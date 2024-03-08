from locidex.classes import run_command

class mafft:
    input_fasta = ''
    params = {}
    alignment = {
        'params':{},
        'alignment':{}
    }
    status = True
    messages = []

    def __init__(self,input_fasta,params):
        (stdout, stderr) = self.align(input_fasta,params)
        self.alignment['params'] = params
        self.alignment['alignment'] = self.parse_align(stdout)
        self.messages.append(stderr)
        if len(self.alignment['alignment']) == 0:
            self.status = False


    def get_alignment(self):
        return self.alignment

    def align(self,input_fasta,params):
        command = ['mafft']
        for p in params:
            command.append(f'--{p}')
            command.append(params[p])
        command.append(input_fasta)
        command = " ".join([str(x) for x in command])
        (stdout, stderr) = run_command(command)
        return (stdout,stderr)


    def parse_align(self,align):
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


class apply_variants:
    seq_record = {
        'id':'',
        'matched_bases':0,
        'num_diff':0,
        'seq':''
    }
    def __init__(self,alignment,ref_id,query_id,fill_gap=False):
        self.call_seq(alignment,ref_id,query_id,fill_gap)

    def trim_ends(self,seq1,seq2):
        seq = []
        is_end = True
        length = len(seq1)
        for i in range(0,length):
            base = seq1[i]
            if is_end and base == '-':
                seq.append('-')
            else:
                is_end = False
                seq.append(seq2[i])

        is_end = True
        length = len(seq1)
        for i in reversed(range(0,length)):
            base = seq1[i]
            if is_end and base == '-':
                seq.append('-')
            else:
                is_end = False
                seq.append(seq2[i])

        return seq

    def extend(self, seq1, seq2):
        pass

    def call_seq(self,alignment,ref_id,query_id,fill_gap=False,):
        ref_seq = alignment[ref_id]
        query_seq = alignment[query_id]
        seq = []
        for idx,base in enumerate(ref_seq):
            if base == '-':
                seq.append('-')
            else:
                if fill_gap and query_seq[idx] == '-':
                    seq.append(base)
                else:
                    seq.append(query_seq[idx])
        seq = "".join([str(x) for x in seq])
        (match, diff ) = self.count_identity(ref_seq,seq)
        self.seq_record['id'] = query_id
        self.seq_record['matched_bases'] = match
        self.seq_record['num_diff'] = diff
        self.seq_record['seq'] = seq

    def count_identity(self,seq1,seq2):
        match = 0
        diff = 0
        for idx, base in enumerate(seq1):
            if seq1[idx] != '-' and seq1[idx] == seq2[idx]:
                match+=1
            else:
                diff+=1
        return [match, diff]





