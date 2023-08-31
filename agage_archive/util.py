import os
import json
import pytz
import yaml

from agage_archive import get_path, open_data_file


def setup():
    """ Setup the config.yml file for the agage_archive package.
    """

    header = '''# Use this file to store configuration settings
# All paths are relative to the network subfolder in the data directory
# If you need to put data files elsewhere, you'll need to use symlinks
---
'''

    config = {}

    # User name
    usr = input("Name (press enter for system ID):") or None
    config["user"] = {"name": usr}

    config["paths"] = {
        "agage_test":
            {
                "agage_md_path": "data-nc",
                "agage_gcms_path": "data-gcms-nc",
                "ale_path": "ale",
                "gage_path": "gage",
                "output_path": "output"
            },
        "agage":
            {
                "agage_md_path": "data-nc.zip",
                "agage_gcms_path": "data-gcms-nc.zip",
                "ale_path": "ale_gage_sio1993/ale",
                "gage_path": "ale_gage_sio1993/gage",
                "output_path": "output"
        }
    }

    with open(get_path('config.yaml'), 'w') as configfile:
        # Write header lines
        configfile.write(header)
        # Dump config dictionary as yaml
        yaml.dump(config, configfile,
                  default_flow_style=False,
                  sort_keys=False)
        
    print(f"Config file written to {get_path('config.yaml')}")
    print("Config file has been populated with default sub-paths relative to data/network. " + \
          "If you want to move the data elsewhere, manually modify the sub-paths in the config file. " + \
          "If the files need to be outside of the data directory, use symlinks.")


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
    
    # Take username from config file if it exists, otherwise try to get it from the system
    with open(get_path("config.yaml"), "r") as f:
        config = yaml.safe_load(f)

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
