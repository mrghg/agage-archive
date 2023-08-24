import configparser
from pathlib import Path
import os
import json
import pytz

from agage_archive import get_path
from agage_archive.io import open_data_file

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
    
    # User name
    usr = input("Name (press enter for system ID):") or None
    config["User"] = {"name": usr}

    with open(get_path('config.ini'), 'w') as configfile:
        config.write(configfile)


def is_number(s):
    """ Check if a string is a number. 
    
    Args:
        s (str): String to check

    Returns:
        bool: True if s is a number, False otherwise
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def lookup_username():
    '''Look up username

    Returns:
        str: Username
    '''

    # Check if config file exists
    config_file = get_path("config.ini")
    if not config_file.exists():
        raise FileNotFoundError(
            "Config file not found. Try running util.setup first")
    
    # Take username from config file if it exists, otherwise try to get it from the system
    config = configparser.ConfigParser()
    config.read(config_file)

    if "User" in config:
        return config["User"]["name"]
    else:
        try:
            return os.environ["USER"]
        except:
            try:
                return os.environ["USERNAME"]
            except:
                try:
                    return os.environ["LOGNAME"]
                except:
                    return "unknown user"


def tz_local_to_utc(index, network, site):
    """ Convert local time to UTC. 
    
    Args:
        index (pandas.DatetimeIndex): Datetime index in local time
        network (str): Network name
        site (str): Site name

    Returns:
        pandas.DatetimeIndex: Datetime index in UTC
    """

    with open_data_file("ale_gage_sites.json", network=network) as file:
        site_info = json.load(file)

    tzoffset_hours = site_info[site]["tz"].split("UTC")[1]

    local_offset = pytz.FixedOffset(int(tzoffset_hours)*60)
    
    ind = index.tz_localize(local_offset)
    
    return ind.tz_convert(None)


if __name__ == "__main__":
    setup()
