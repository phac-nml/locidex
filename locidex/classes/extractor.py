import pandas as pd
import numpy as np
from locidex.constants import NT_SUB, STOP_CODONS, START_CODONS
class extractor:
    seqs = {}
    df = pd.DataFrame()
    def __init__(self,df,seq_data,sseqid_col,queryid_col,qstart_col,qend_col,qlen_col,sstart_col,send_col,slen_col,sstrand_col,bitscore_col,overlap_thresh=1,extend_threshold_ratio = 0.2,filter_contig_breaks=True):
        self.filter_contig_breaks = filter_contig_breaks
        self.df = self.set_extraction_pos(df, sstart_col, send_col)
        self.is_complete(self.df,qstart_col,qend_col,qlen_col)
        self.is_contig_boundary(self.df,'ext_start','ext_end',slen_col)
        if filter_contig_breaks:
            self.df = df[ (df['is_5prime_boundary'] == False) & (df['is_3prime_boundary'] == False)]
            self.df = self.df.reset_index(drop=True)
        self.set_revcomp(self.df,sstart_col,send_col,sstrand_col)
        pcols = [qstart_col,qend_col,sstart_col,send_col]
        for c in pcols:
            self.df[c] = self.df[c].apply(lambda x: x - 1)
        #self.df = self.fix_postioning(self.df)

        sort_cols = [sseqid_col, queryid_col, sstart_col, send_col, bitscore_col]
        ascending_cols = [True, True, True, True, False]
        #self.df = self.recursive_filter_overlap_records(self.df, sseqid_col, bitscore_col, sort_cols, ascending_cols, overlap_threshold=overlap_thresh)
        self.df = self.extend(self.df,sseqid_col, queryid_col, qstart_col, qend_col, sstart_col,send_col,slen_col, qlen_col, bitscore_col, overlap_threshold=1)

        #self.df = self.recursive_filter_overlap_records(self.df, sseqid_col, bitscore_col, sort_cols, ascending_cols,
        #                                                overlap_threshold=overlap_thresh)
        self.df = self.set_extraction_pos(self.df, sstart_col, send_col)
        loci_ranges = self.group_by_locus(self.df,sseqid_col, queryid_col,qlen_col,extend_threshold_ratio)
        self.seqs = self.extract_seq(loci_ranges, seq_data)
        pass

    def is_5prime_complete(self,df,qstart_col):
        df['is_5prime_complete'] = np.where(df[qstart_col] == 1, True, False)

    def is_3prime_complete(self,df,qend_col,qlen_col):
        df['is_3prime_complete'] = np.where(df[qend_col] == df[qlen_col] , True, False)

    def is_complete(self,df,qstart_col,qend_col,qlen_col):
        self.is_5prime_complete(df,qstart_col)
        self.is_3prime_complete(df,qend_col,qlen_col)
        df['is_complete'] = np.where(((df['is_5prime_complete'] == True) &  (df['is_3prime_complete'] == True)), True, False)

    def is_contig_boundary(self,df,sstart_col,send_col,slen_col):
        df['is_5prime_boundary'] = np.where(df[sstart_col] == 1, True, False)
        df['is_3prime_boundary'] = np.where(df[send_col] == df[slen_col], True, False)
        df['is_on_boundary'] = np.where(((df['is_5prime_boundary'] == True) & (df['is_3prime_boundary'] == True)), True, False)

    def set_revcomp(self,df,sstart_col,send_col,strand_col):
        df['reverse'] = np.where(df[sstart_col] > df[send_col], True, False)
        df['complement'] = np.where(df[strand_col] == 'plus', True, False)

    def set_extraction_pos(self,df,start_col,end_col):
        df['ext_start'] = df[start_col]
        df['ext_end'] = df[end_col]
        for idx, row in df.iterrows():
            start = row['ext_start']
            end = row['ext_end']
            if start > end:
                tmp = end
                end = start
                start = tmp
                df.loc[idx, 'ext_start'] = start
                df.loc[idx, 'ext_end'] = end
        return df

    def fix_postioning(self,df):
        for idx, row in df.iterrows():
            start = row['ext_start']
            end = row['ext_end']
            rev = row['reverse']
            complement = row['complement']

            if rev and complement:
                start += 1
                end += 1
            elif rev:
                start+=1
            elif complement:
                end += 1
            else:
                end += 1

            if row['reverse'] == True:
                tmp = row['ext_end']
                end = row['ext_start']
                start = tmp
            df.loc[idx, 'ext_start'] = start
            df.loc[idx, 'ext_end'] = end

        return df

    def remove_redundant_hits(self,df,seqid_col, bitscore_col, overlap_threshold=1):
        seq_id_list = list(df[seqid_col].unique())
        filter_df = []
        for seqid in seq_id_list:
            subset = df[df[seqid_col] == seqid]
            prev_contig_id = ''
            prev_index = -1
            prev_contig_start = -1
            prev_contig_end = -1
            prev_score = 0
            filter_rows = []
            for idx, row in subset.iterrows():
                contig_id = row[seqid_col]
                contig_start = row['ext_start']
                contig_end = row['ext_end']
                score = float(row[bitscore_col])

                if prev_contig_id == '':
                    prev_index = idx
                    prev_contig_id = contig_id
                    prev_contig_start = contig_start
                    prev_contig_end = contig_end
                    prev_score = score
                    continue

                if (contig_start >= prev_contig_start and contig_start <= prev_contig_end) or (
                        contig_end >= prev_contig_start and contig_end <= prev_contig_end):
                    overlap = abs(contig_start - prev_contig_end)

                    if overlap > overlap_threshold:
                        if prev_score < score:
                            filter_rows.append(prev_index)
                        else:
                            filter_rows.append(idx)

                prev_index = idx
                prev_contig_id = contig_id
                prev_contig_start = contig_start
                prev_contig_end = contig_end
                prev_score = score


            valid_ids = list( set(subset.index) - set(filter_rows)  )

            filter_df.append(subset.filter(valid_ids, axis=0))


        return pd.concat(filter_df, ignore_index=True)


    def recursive_filter_overlap_records(self,df, seqid_col, bitscore_col, sort_cols, ascending_cols, overlap_threshold=1):
        size = len(df)
        prev_size = 0
        while size != prev_size:
            df = df.sort_values(sort_cols,
                                ascending=ascending_cols).reset_index(drop=True)
            df = self.remove_redundant_hits(df, seqid_col, bitscore_col, overlap_threshold=overlap_threshold)
            prev_size = size
            size = len(df)

        return df.sort_values(sort_cols,ascending=ascending_cols).reset_index(drop=True)

    def extract_seq(self,loci_data,seq_data):
        seqs = []
        id = 0
        for locus_name in loci_data:
            for row in loci_data[locus_name]:
                query_id = row['query_id']
                start = row['start']
                end = row['end'] + 1
                seqid = row['seqid']
                is_reverse = row['reverse']
                is_complement = row['complement']
                is_extended = row['is_extended']


                is_complete = row['is_complete']
                fivep_trunc = row['is_5prime_boundary']
                threep_trunc = row['is_3prime_boundary']
                is_trunc = False
                if fivep_trunc or threep_trunc:
                    is_trunc = True

                if seqid in seq_data:

                    seq = seq_data[seqid]['seq'][start:end]
                    if is_reverse and not is_complement:
                        seq = seq[::-1].translate(NT_SUB)
                    start_codon = seq[0:3]
                    stop_codon = seq[-3:]
                    if stop_codon in STOP_CODONS:
                        is_stop_valid = True
                    else:
                        is_stop_valid = False

                    if start_codon in START_CODONS:
                        is_start_valid = True
                    else:
                        is_start_valid = False

                    if is_start_valid and is_stop_valid:
                        cds_valid = True
                    else:
                        cds_valid = False


                    seqs.append({'seqid':seqid, 'id':str(id),
                                 'locus_name':locus_name,'query_id':query_id,
                                 'start':start, 'end':end,
                                 'reverse':is_reverse,
                                 'complement':is_complement,
                                 'is_complete':is_complete,
                                 'is_trunc':is_trunc,
                                 'fivep_trunc':fivep_trunc,
                                 'threep_trunc':threep_trunc,
                                 'seq':seq,'start_codon':start_codon,'stop_codon':stop_codon,
                                 'is_stop_valid':is_stop_valid,'is_start_valid':is_start_valid,
                                 'is_cds_valid':cds_valid})
                    id+=1
        return seqs

    def extend(self,df,seqid_col, queryid_col, qstart_col, qend_col, sstart_col,send_col,slen_col, qlen_col, bitscore_col, overlap_threshold=1):
        sort_cols = [seqid_col,'ext_start', 'ext_end', bitscore_col]
        ascending_cols = [True, True, True, False]
        df = self.recursive_filter_overlap_records(df, seqid_col, bitscore_col, sort_cols, ascending_cols, overlap_threshold)
        df = df.sort_values(sort_cols,
                            ascending=ascending_cols).reset_index(drop=True)

        queries = df[queryid_col].to_list()

        #Remove incomplete hits when complete ones are present
        filtered = []
        for query in queries:
            subset = df[df[queryid_col] == query]
            complete = subset[subset['is_complete'] == True]
            num_complete = len(complete)
            if num_complete > 0:
                filtered.append(complete)
            else:
                filtered.append(subset)
        df = pd.concat(filtered, ignore_index=True)

        trunc_records = df[df['is_complete'] == False]
        if len(trunc_records) == 0:
            return df

        is_extended = []
        five_p_ext = []
        three_p_ext = []
        #Extend truncated sequences
        for idx, row in df.iterrows():
            e = False
            fivep_e = False
            threep_e = False
            if row['is_complete']:
                is_extended.append(e)
                five_p_ext.append(fivep_e)
                three_p_ext.append(threep_e)
                continue

            qstart = int(row[qstart_col])
            qend = int(row[qend_col])
            qlen = int(row[qlen_col])
            sstart = int(row[sstart_col])
            send = int(row[send_col])
            slen = int(row[slen_col])
            is_rev = row['reverse']

            five_p_complete = row['is_5prime_complete']
            five_p_delta = qstart
            three_p_complete = row['is_3prime_complete']
            three_p_delta = qlen - qend
            if three_p_delta < 0:
                three_p_delta = 0


            if not five_p_complete:
                e = True
                fivep_e = True
                if is_rev:
                    sstart += five_p_delta
                else:
                    sstart -= five_p_delta

            if sstart < 1:
                sstart = 0

            if not three_p_complete:
                e = True
                threep_e = True
                if is_rev:
                    send -= three_p_delta
                else:
                    send += three_p_delta

            if send >= slen:
                send = slen - 1

            is_extended.append(e)
            five_p_ext.append(fivep_e)
            three_p_ext.append(threep_e)


            row[qstart_col] = qstart
            row[qend_col] = qend
            row[qlen_col] = qlen
            row[sstart_col] = sstart
            row[send_col] = send
            row[slen_col] = slen

            df.loc[idx] = row

        df['is_extended'] = is_extended
        df['is_5p_extended'] = five_p_ext
        df['is_3p_extended'] = three_p_ext
        return df

    def group_by_locus(self,df,seqid_col,query_col,qlen_col,extend_threshold_ratio = 0.2):
        sort_cols = ['locus_name',query_col,'ext_start', 'ext_end']
        ascending_cols = [True, True, True, True]
        df = df.sort_values(sort_cols,
                            ascending=ascending_cols).reset_index(drop=True)

        queries = df['qseqid'].to_list()

        #Remove incomplete hits when complete ones are present
        filtered = []
        for query in queries:
            subset = df[df['qseqid'] == query]
            complete = subset[subset['is_complete'] == True]
            num_complete = len(complete)
            if num_complete > 0:
                filtered.append(complete)
            else:
                filtered.append(subset)
        df = pd.concat(filtered, ignore_index=True)

        sort_cols = ['locus_name',seqid_col,'ext_start', 'ext_end']
        ascending_cols = [True, True, True, True]
        df = df.sort_values(sort_cols,
                            ascending=ascending_cols).reset_index(drop=True)

        loci = {}
        for idx, row in df.iterrows():
            query_id = row[query_col]
            locus_name = row['locus_name']
            start = row['ext_start']
            end = row['ext_end']
            qlen = row[qlen_col]
            seqid = row[seqid_col]
            is_extended = row['is_extended']
            is_reverse = row['reverse']
            is_complement = row['complement']
            is_complete = row['is_complete']
            is_5prime_boundary = row['is_5prime_boundary']
            is_3prime_boundary = row['is_3prime_boundary']
            is_5p_extended = row['is_5p_extended']
            is_3p_extended = row['is_3p_extended']
            #if not is_complement and is_reverse:
            #    start-=1

            if locus_name in loci:
                found = False
                for i in range(0,len(loci[locus_name])):
                    s = loci[locus_name][i]['start']
                    e = loci[locus_name][i]['end']
                    sid = loci[locus_name][i]['seqid']
                    if sid != seqid:
                        continue
                    if (s >= start and s <= end) or (e >= start and e<= end):
                        found = True
                        if end > e:
                            loci[locus_name][i]['end'] = end
                            break
                    if s >= end:
                        gap = s - end
                        max_gap = qlen * extend_threshold_ratio
                        if gap <= max_gap:
                            loci[locus_name][i]['end'] = end
                            found = True
                            break
                if found:
                    continue
            else:
                loci[locus_name] = []

            loci[locus_name].append({
                'seqid':seqid,
                'query_id':query_id,
                'start':start,
                'end':end,
                'reverse':is_reverse,
                'complement':is_complement,
                'is_complete':is_complete,
                'is_extended':is_extended,
                'is_5p_extended':is_5p_extended,
                'is_3p_extended': is_3p_extended,
                'is_5prime_boundary':is_5prime_boundary,
                'is_3prime_boundary':is_3prime_boundary

            })

        return loci



