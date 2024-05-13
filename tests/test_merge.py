import pytest
from locidex import merge


DUPLICATE_NAMES = [
    "locidex/example/merge/merge_in/report.json", 
    "locidex/example/merge/merge_in/report2.json", 
]

MERGE_SUCCESSFULLY = [
    "locidex/example/merge/merge_in/report.json", 
    "locidex/example/merge/merge_in/report1.json", 
]

def test_throws_duplicate_error():
    """
    """
    with pytest.raises(SystemExit, match="Attempting to merge allele profiles with the same sample name: 1"):
        merge.read_file_list(DUPLICATE_NAMES, perform_validation=True)

def test_merge_profiles_no_error():
    merge.read_file_list(MERGE_SUCCESSFULLY, perform_validation=True)
    