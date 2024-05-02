"""
Test the ever growing util functions
"""

import pytest
from locidex import utils
from locidex import manifest
from argparse import Namespace


def test_check_db_groups_pass(monkeypatch):
    nm_group = Namespace(db_group="Db1", db_name="test_name", db_version="1.0.0")
    analysis_params = {"db_group": "Db1", "db_name": "test_name", "db_version": "1.0.0"}

    def mockreturn(*args, **kwargs):
        return True
    monkeypatch.setattr(manifest, "get_manifest_db", mockreturn)
    analysis_params = utils.check_db_groups(analysis_params, nm_group)
    assert analysis_params["db"]

def test_check_db_groups_fail():
    nm_group = Namespace(db_group="Db1", db_name="test_name", db_version="1.0.0")
    analysis_params = {"db_group": "Db1", "db_name": "test_name"}

    with pytest.raises(KeyError):
        analysis_params = utils.check_db_groups(analysis_params, nm_group)