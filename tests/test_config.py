from pathlib import Path
import yaml

from agage_archive.config import Paths


repo_path = Path(__file__).resolve().parents[1]


def test_config_exists():
    """ Check that the configuration file exists
    This is just so to help diagnose this as a problem if subsequent tests fail
    """

    config_file = repo_path / "agage_archive" / "config.yaml"

    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found at {config_file}. Create one before running tests")

    # Read the configuration file
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Check that agage-test network is in the configuration file
    if not "agage_test" in config["paths"].keys():
        raise KeyError("agage_test network not found in configuration file. Must exist to run tests")


def test_paths():

    # Check that "this_repo" option returns the repository path
    path_this_repo = Paths(this_repo=True)
    assert path_this_repo.root == repo_path / "agage_archive"
    assert path_this_repo.data == repo_path / "data"

    # Read the configuration file for agage_test network
    config_file = repo_path / "agage_archive" / "config.yaml"
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Check that some of the agage_test paths are as expected
    path = Paths("agage_test")
    assert path.output_path == \
              config["paths"]["agage_test"]["output_path"]

    # If we try to retrieve a path that doesn't exist, it should raise a KeyError
    try:
        path = Paths("non_existent_network")
    except KeyError:
        pass
    else:
        raise AssertionError("KeyError not raised for non-existent network")

    # If we try to retrieve a path that doesn't exist with errors ignored, it should not error
    path = Paths("non_existent_network", errors="ignore")

