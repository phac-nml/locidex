import pytest
import pandas as pd
from locidex.classes.assignment import assignment

# Fixture to create a sample DataFrame
@pytest.fixture
def sample_dataframe():
    data = {
        'query_id': [1, 2, 2, 3],
        'value': [100, 200, 200, 300],
        'another_column': ['A', 'B', 'B', 'C']
    }
    df = pd.DataFrame(data)
    loci_metadata = {'locus_1': 'metadata_1', 'locus_2': 'metadata_2'}
    return df, loci_metadata

def test_assignment_initialization(sample_dataframe):
    df, loci_metadata = sample_dataframe
    query_id_col = 'query_id'
    assgn = assignment(df, query_id_col, loci_metadata)
    
    # Check if the DataFrame is correctly set
    assert assgn.df.equals(df), "DataFrame not correctly assigned."
    
    # Check if the columns attribute correctly lists the DataFrame's columns
    expected_columns = ['query_id', 'value', 'another_column']
    assert assgn.columns == expected_columns, "Columns attribute not correctly set."
    
    # Check if loci_metadata is correctly set
    expected_loci_metadata = {'locus_1': 'metadata_1', 'locus_2': 'metadata_2'}
    assert assgn.loci_metadata == expected_loci_metadata, "Loci metadata not correctly assigned."
    
    # Check if the unique queries list is correctly set
    expected_queries = [1, 2, 3]
    assert assgn.queries == expected_queries, "Queries list not correctly set."