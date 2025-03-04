import pytest
from locidex import merge
import errno
import json


DUPLICATE_NAMES = [
    "locidex/example/merge/merge_in/report.json",
    "locidex/example/merge/merge_in/report2.json",
]

MERGE_SUCCESSFULLY = [
    "locidex/example/merge/merge_in/report.json",
    "locidex/example/merge/merge_in/report1.json",
]

MERGE_MISMATCH_PROFILES = "locidex/example/merge/merge_inputassure/sample1.mlst.json"

def test_throws_duplicate_error():
    """
    """
    with pytest.raises(ValueError, match="Attempting to merge allele profiles with the same sample name: 1"):
        merge.read_file_list(DUPLICATE_NAMES, perform_db_validation=True)

def test_merge_profiles_no_error():
    merge.read_file_list(MERGE_SUCCESSFULLY, perform_db_validation=True)


def test_check_files_exist():
    with pytest.raises(FileNotFoundError, match=str(errno.ENOENT)):
        fail = "dne.txt"
        merge.check_files_exist([fail])

    merge.check_files_exist(MERGE_SUCCESSFULLY)

def test_extract_profiles():
    records, compare_error = merge.read_file_list(MERGE_SUCCESSFULLY, perform_db_validation=True)
    extracted_profiles = merge.extract_profiles(records)
    assert len(extracted_profiles) == 2
    value1, value2 = extracted_profiles.values()
    assert value1 == value2
    assert value1 == {'aroC': '9048803cd72dee3c868cd2dc5dc5650d', 'dnaN': '2772ad8b8e0f7b50f1396c31fbe53f2d', 'hemD': '620f99723c4e190abe096b11ca34b944', 'hisD': '38027ac1ac34817584a176c7e575e97e', 'purE': '9855cbf4009439498bf84cacefce4d8f', 'sucA': '9289fc07cc8e93cfe0716e6f613cefdb', 'thrA': '9e1aa76bb42279ed7ec8fc30f984b65d'}

def test_input_assure():
    mlst_file = json.load(open(MERGE_MISMATCH_PROFILES))
    ## Test that if the sample name matches the profile in the MLST json file nothing is changed
    new_mlst, mlst_report = merge.validate_profiles(mlst_file, "sampleA", "sample1.mlst.json")
    assert new_mlst["data"]["profile"] == {'sampleA': {'l1': '1', 'l2': '1', 'l3': '1'}}
    assert mlst_report == None
    ## Test that a different sample than the profile in the MLST json will be changed
    same_mlst, mlst_report = merge.validate_profiles(mlst_file, "sample1", "sample1.mlst.json")
    assert same_mlst["data"]["profile"] == {'sample1': {'l1': '1', 'l2': '1', 'l3': '1'}}
    assert mlst_report[2] == "sample1 ID and JSON key in sample1.mlst.json DO NOT MATCH. The 'sampleA' key in sample1.mlst.json has been forcefully changed to 'sample1': User should manually check input files to ensure correctness."




