from Bio import GenBank
import gzip
from mimetypes import guess_type
from functools import partial
import os
from locidex.utils import revcomp,calc_md5

class parse_gbk:
    input_file = None
    seq_obj = None
    status = True
    messages = []

    def __init__(self,input_file):
        self.input_file= input_file
        if not os.path.isfile(self.input_file):
            self.messages.append(f'Error {self.input_file} does not exist')
            self.status = False
        self.seq_obj = self.parse_reference_gbk()


    def get_acs(self):
        if self.seq_obj is not None:
            return list(self.seq_obj.keys())
        else:
            return []

    def get_seq_by_acs(self,acs):
        if self.seq_obj is not None and acs in self.seq_obj:
            return self.seq_obj[acs]
        else:
            return {}

    def get_feature(self,acs,feature):
        s = self.get_seq_by_acs(acs)
        if feature in s['features']:
            return s['features'][feature]
        else:
            return {}


    def parse_reference_gbk(self):
        """
        Method to parse the GenBank reference file, clean the strings, and
        return the reference features of interest.

        Parameters
        ----------
        gbk_file : str
            The file path to the reference genbank format file with sequence annotations
        Returns
        -------
        dict
            A dictionary of all the reference features
        """
        encoding = guess_type(self.input_file)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        sequences = {}
        with _open(self.input_file) as handle:
            for record in GenBank.parse(handle):
                gb_accession = record.locus
                gb_accession_version = 1
                if len(record.accession) == 2:
                    gb_accession = record.accession[0]
                    gb_accession_version = gb_accession[1]
                # clean the sequence
                genome_seq = repr(record.sequence).replace("\'", '')
                sequences[gb_accession] = {
                    'accession': gb_accession,
                    'version': gb_accession_version,
                    'features': {'source': genome_seq}
                }
                features = record.features
                # retrieve the features if present
                for feat in features:
                    if feat.key == 'CDS' or feat.key == '5\'UTR' or feat.key == '3\'UTR':
                        if not feat.key in sequences[gb_accession]['features']:
                            sequences[gb_accession]['features'][feat.key] = []
                        qualifier = feat.qualifiers
                        positions = []
                        gene_name = ''
                        locus_tag = ''
                        aa = ''
                        # more string cleaning
                        for name in qualifier:
                            if name.key == '/gene=':
                                gene_name = name.value.replace("\"", '').strip()
                            if name.key == '/translation=':
                                aa = name.value.replace("\"", '').strip()
                            if name.key == '/locus_tag=':
                                gene_name = name.value.replace("\"", '').strip()
                                locus_tag = gene_name
                        if locus_tag != '':
                            gene_name = locus_tag
                        locations = feat.location.strip().replace("join(", '').replace(')', '').split(',')
                        seq = []
                        # retreive the locations of the features
                        for location in locations:
                            # more string cleaning
                            location = location.replace('<', '').replace('>', '')
                            if not 'complement' in location:
                                location = location.split('.')
                                start = int(location[0]) - 1
                                end = int(location[2])
                                seq.append(genome_seq[start:end].replace("\'", ''))
                                positions.append([start, end])
                            else:
                                location = location.replace('complement(', '').replace(')', '').split('.')
                                start = int(location[0]) - 1
                                end = int(location[2])
                                seq.append(revcomp(genome_seq[start:end].replace("\'", '')))
                                positions.append([start, end])

                        seq = ''.join(seq)
                        sequences[gb_accession]['features'][feat.key].append(
                            {'gene_name': gene_name, 'dna_seq': seq.lower(), 'aa_seq': aa.lower(), 'positions': positions,
                             'gene_len': len(seq),'aa_hash':calc_md5([aa])[0],'dna_hash':calc_md5([seq])[0]})

        return sequences