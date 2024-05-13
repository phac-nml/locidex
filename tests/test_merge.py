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
    value1, value2 = extracted_profiles.values()
    assert value1 == value2
    assert value1 == {'aroC': '9048803cd72dee3c868cd2dc5dc5650d', 'dnaN': '2772ad8b8e0f7b50f1396c31fbe53f2d', 'hemD': '620f99723c4e190abe096b11ca34b944', 'hisD': '38027ac1ac34817584a176c7e575e97e', 'purE': '9855cbf4009439498bf84cacefce4d8f', 'sucA': '9289fc07cc8e93cfe0716e6f613cefdb', 'thrA': '9e1aa76bb42279ed7ec8fc30f984b65d'}