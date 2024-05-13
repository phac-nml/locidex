import pytest
from locidex import merge
import errno


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
    

def test_check_files_exist():
    with pytest.raises(SystemExit, match=str(errno.ENOENT)):
        fail = "dne.txt"
        merge.check_files_exist([fail])
    
    merge.check_files_exist(MERGE_SUCCESSFULLY)

def test_extract_profiles():
    records = merge.read_file_list(MERGE_SUCCESSFULLY, perform_validation=True)
    extracted_profiles = merge.extract_profiles(records)
    assert len(extracted_profiles) == 2
    key1, key2 = extracted_profiles.keys()
    assert extracted_profiles[key1] == extracted_profiles[key2]