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

    def __init__(self,input_fasta,params={}):
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






