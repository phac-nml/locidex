import json
import pytest
import tempfile
import os
from pathlib import Path
from locidex.classes.db import db_config, search_db_conf
import locidex.constants as constants 

required_fields = ["database", "user", "password"]

# Fixture for creating and cleaning up a temporary configuration file
@pytest.fixture
def temp_config():
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w')
    yield temp_file.name  # Provide the path to the test
    os.unlink(temp_file.name)  # Cleanup

def write_temp_config(path, data):
    with open(path, 'w') as tmp:
        json.dump(data, tmp)

def test_db_config_valid(temp_config):
    valid_config_data = {
        "database": "test_db",
        "user": "test_user",
        "password": "test_pass"
    }
    write_temp_config(temp_config, valid_config_data)
    config = db_config(temp_config, required_fields)
    assert config.status is True
    assert config.messages == []

def test_db_config_empty_file(temp_config):
    # Assuming that the locidex errors out with Json decode error if JSON file is empty:
    with pytest.raises(json.JSONDecodeError):
        db_config(temp_config, required_fields)

    # config = db_config(temp_config, required_fields)
    # assert config.status is False
    # assert "required field missing" in config.messages[0]

def test_db_config_no_required_fields_present(temp_config):
    no_required_fields_data = {"unrelated_field": "some_value"}
    write_temp_config(temp_config, no_required_fields_data)
    config = db_config(temp_config, required_fields)
    assert config.status is False
    assert "required field missing" in config.messages[0]

def test_db_config_initialization(temp_config):
    # Create a configuration with all required fields present
    config_data = {"database": "test_db", "user": "test_user", "password": "test_pass"}
    write_temp_config(temp_config, config_data)
    
    config = db_config(temp_config, required_fields)
    
    assert config.status is True, "Initialization failed when it should have succeeded."
    assert config.config == config_data, "Configuration data was not correctly parsed."

def test_validate_dict_keys():
    # Create a minimal but valid configuration file
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        # Assuming an empty JSON object is a valid configuration for your class
        json.dump({}, tmp)
        tmp_path = tmp.name

    # Now that we have a valid (though minimal) configuration file, proceed with testing
    try:
        config = db_config(tmp_path, required_fields=[])

        # Test validation logic
        data_dict = {"database": "db", "user": "usr", "password": "pwd"}
        status, missing_key = config.validate_dict_keys(data_dict, ["database", "user", "password"])
        assert status is True, "Failed to validate correctly when all keys are present."
        assert missing_key == '', "Reported a missing key when none was missing."

        status, missing_key = config.validate_dict_keys(data_dict, ["missing_key"])
        assert status is False, "Failed to identify missing key."
        assert missing_key == "missing_key", "Identified the wrong missing key."
    finally:
        # Cleanup - remove the temporary file
        os.remove(tmp_path)

## search db config tests

@pytest.fixture
def temp_db_dir(tmp_path):
    # Create a temporary directory structure for testing
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    config_file = db_dir / "db_config.json"
    meta_file = db_dir / "db_meta.json"
    blast_dir = db_dir / "blast"
    blast_dir.mkdir()
    nucleotide_dir = blast_dir / "nucleotide"
    nucleotide_dir.mkdir()
    protein_dir = blast_dir / "protein"
    protein_dir.mkdir()
    
    # Create mock files
    config_file.write_text(json.dumps({i: "value" for i in constants.DBConfig._keys()}))
    meta_file.write_text(json.dumps({"meta": "data"}))
    (nucleotide_dir / "nucleotide").touch()
    (protein_dir / "protein").touch()
    
    yield db_dir

def test_search_db_conf_initialization_and_blast_paths_setup(temp_db_dir):
    # Given that the intialization of the class performs setup operations, the init and blast path setup were tested together
    # Setup basenames and required fields
    db_basenames = {
        "config": "db_config.json",
        "meta": "db_meta.json"
    }
    #required_fields = ["key"]
    #required_fields = [*constants.DB_CONFIG_FIELDS]
    required_fields = [*constants.DBConfig._keys()]
    
    # Initialize search_db_conf
    search_conf = search_db_conf(str(temp_db_dir), db_basenames, required_fields)
    
    # Assertions to verify initialization was successful
    assert search_conf.status is True, "Failed to initialize with valid directory and files"
    assert search_conf.config_file_path == str(temp_db_dir / "db_config.json"), "Config file path incorrect"
    assert search_conf.meta_file_path == str(temp_db_dir / "db_meta.json"), "Meta file path incorrect"

    # Further assertions to verify blast paths setup as part of the initialization
    assert search_conf.blast_paths['nucleotide'] == str(temp_db_dir / "blast/nucleotide/nucleotide"), "Nucleotide BLAST path incorrect"
    assert search_conf.blast_paths['protein'] == str(temp_db_dir / "blast/protein/protein"), "Protein BLAST path incorrect"
