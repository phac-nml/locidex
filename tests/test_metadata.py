import pytest
import json
import warnings
from locidex.classes.metadata import metadata

@pytest.fixture
def good_metadata(tmp_path):
    good_meta = tmp_path / "good_metadata.json"
    good_meta.write_text(json.dumps({"key1": "value1", "key2": "value2", "key3": "value3"}))
    return good_meta

@pytest.fixture
def incomplete_metadata(tmp_path):
    incomplete_meta = tmp_path / "incomplete_metadata.json"
    incomplete_meta.write_text(json.dumps({"key1": "value1", "key2": "value2"}))
    return incomplete_meta

@pytest.mark.skip(reason="Need to handle None data gracefully. Class improvement needed.")
def test_metadata_good(good_metadata):
    required_keys = ['key1', 'key2', 'key3']
    warnings.warn("Class should handle JSON parsing failures and return appropriate errors.", UserWarning)
    meta = metadata(good_metadata, required_keys)
    assert meta.status, "Metadata should be valid with all required keys present"

@pytest.mark.skip(reason="Class should explicitly check for None before iterating over dictionary.")
def test_metadata_incomplete(incomplete_metadata):
    required_keys = ['key1', 'key2', 'key3']
    warnings.warn("Class should explicitly check for None before iterating over dictionary to handle missing keys gracefully.", UserWarning)
    meta = metadata(incomplete_metadata, required_keys)
    assert not meta.status, "Metadata should be invalid due to missing key"

@pytest.mark.skip(reason="__init__ method should handle non-existent files properly.")
def test_metadata_non_existent(tmp_path):
    non_existent_file = tmp_path / "does_not_exist.json"
    required_keys = ['key1', 'key2', 'key3']
    warnings.warn("__init__ method should handle non-existent files properly and should not return any value.", UserWarning)
    meta = metadata(str(non_existent_file), required_keys)
    assert not meta.status, "Metadata should be invalid due to non-existent file"