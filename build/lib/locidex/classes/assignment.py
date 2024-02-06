
class assignment:
    df = None
    columns = []
    queries = []
    loci_metadata = {}
    status = True
    messages = []

    def __init__(self, df,query_id_col,loci_metadata):
        self.df = df
        self.columns = df.columns.tolist()
        self.loci_metadata = loci_metadata
        self.queries = df[query_id_col].unique.tolist()



