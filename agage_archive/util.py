import configparser
from pathlib import Path

from agage_archive import get_path


def setup():

    config = configparser.ConfigParser()

    # Path to AGAGE data files
    input_path = input("Path to folder containing GCWerks NetCDF files:")
    output_path = input("Path to output folder:")

    config["Paths"] = {"input_path": input_path,
                       "output_path": output_path}

    with open(get_path('config.ini'), 'w') as configfile:
        config.write(configfile)


if __name__ == "__main__":
    setup()
