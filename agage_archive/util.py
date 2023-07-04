import configparser
from pathlib import Path

from agage_archive import get_path


def setup():
    """ Setup the config file for the agage_archive package. """

    config = configparser.ConfigParser()

    # Get path to ALE/GAGE data included in repo
    repo_path = Path(__file__).parents[1].absolute()
    ale_default = repo_path / "data/ale_gage_sio1993/ale"
    gage_default = repo_path / "data/ale_gage_sio1993/gage"


    # Paths to data files
    agage_path = input("Path to folder containing AGAGE GCWerks NetCDF files:")
    ale_path = input(f"Path to folder containing ALE tar.gz files " + \
                     f"(press enter for default {ale_default}): ") or ale_default
    gage_path = input("Path to folder containing GAGE tar.gz files:"  + \
                     f"(press enter for default {gage_default}): ") or gage_default
    output_path = input("Path to output folder:")

    config["Paths"] = {"agage_path": agage_path,
                       "ale_path": ale_path,
                       "gage_path": gage_path,
                       "output_path": output_path}

    with open(get_path('config.ini'), 'w') as configfile:
        config.write(configfile)


def is_number(s):
    """ Check if a string is a number. """
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    setup()
