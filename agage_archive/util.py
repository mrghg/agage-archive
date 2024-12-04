import os
import json
import pytz
import yaml
import pandas as pd

from agage_archive.config import Paths, open_data_file, data_file_path


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
    
    paths = Paths()

    # Take username from config file if it exists, otherwise try to get it from the system
    with open(paths.root / "config.yaml", "r") as f:
        config = yaml.safe_load(f)

    if config["user"]["name"] != "":
        return config["user"]["name"]
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


def excel_to_csv(file, network):
    """ Convert Excel data specification file to CSVs.

    Args:
        file (str): File name. Must be either:
            "data_release_schedule", "data_combination" or "data_exclude"
        network (str): Network name
    """

    filename = data_file_path(file + ".xlsx", network)

    if not filename.exists():
        raise ValueError(f"Check filename: {filename}")

    csv_folder_name = filename.parent / filename.name.split(".")[0]

    if not csv_folder_name.exists():
        raise ValueError(f"Create folder {csv_folder_name}")

    # Read Excel file and output worksheets to CSVs
    sheet_names = pd.ExcelFile(filename).sheet_names

    for sheet in sheet_names:
        # Read header
        header_sheet = pd.read_excel(filename,
                                    sheet_name=sheet)

        header = [header_sheet.columns[0]]
        i=0
        while header_sheet.iloc[i, 0][0] == "#":
            header.append(header_sheet.iloc[i, 0])
            i+=1

        # Special case of general release date
        if header_sheet.iloc[i, 0][:7].upper() == "GENERAL":
            header.append("# " + header_sheet.iloc[i, 0] + ": " + header_sheet.iloc[i, 1])
            i+=1

        xlsx_data = pd.read_excel(filename, sheet_name=sheet,
                                skiprows=i+1, dtype=str)
        
        csv_filename = csv_folder_name / f"{file}_{sheet}.csv"

        with open(csv_filename, "w") as f:
            for h in header:
                f.writelines(h + "\n")
            xlsx_data.to_csv(f, index = None)
