import configparser

from agage_archive import get_path


def setup():

    config = configparser.ConfigParser()

    # Paths to data files
    agage_path = input("Path to folder containing AGAGE GCWerks NetCDF files:")
    ale_path = input("Path to folder containing ALE tar.gz files:")
    gage_path = input("Path to folder containing GAGE tar.gz files:")
    output_path = input("Path to output folder:")

    config["Paths"] = {"agage_path": agage_path,
                       "ale_path": ale_path,
                       "gage_path": gage_path,
                       "output_path": output_path}

    with open(get_path('config.ini'), 'w') as configfile:
        config.write(configfile)


if __name__ == "__main__":
    setup()
