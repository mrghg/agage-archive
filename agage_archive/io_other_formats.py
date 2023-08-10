import pandas as pd
from pathlib import Path
import numpy as np
import json

from agage_archive import Paths as Pth


paths = Pth()
home = Path.home()

species_wang = {
    "cfc-11": "CFC-11S",
    "cfc-12": "CFC-12",
    "cfc-113": "CFC-113",
    "ccl4": "CCl4",
    "ch3ccl3": "CH3CCl3",
    "n2o": "N2O",
    "ch4": "CH4"
}

sites_wang = {
    "CMO": "ORG",
    "CGO": "CGO",
    "MHD": "MHD",
    "SMO": "SMO",
    "THD": "THD",
    "ADR": "ADR",
    "RPB": "RPB",
}

def read_wang_file(file):

    # Reading the first 10 lines of the file to understand its structure
    with open(file, "r") as f:
        first_lines = [f.readline().strip() for _ in range(10)]

    # 1. Read the file into a Pandas dataframe with correct headers
    complete_data = pd.read_csv(file, skiprows=6, delim_whitespace=True,
                                names=first_lines[5].split())

    # 2. Construct the datetime column and set it as the index
    complete_data['datetime'] = pd.to_datetime(complete_data['YYYY'].astype(str) + '-' + 
                                            complete_data['MM'].astype(str).str.zfill(2) + '-' + 
                                            complete_data['DD'].astype(str).str.zfill(2) + ' ' + 
                                            complete_data['hh'].astype(str).str.zfill(2) + ':' + 
                                            complete_data['min'].astype(str).str.zfill(2))
    complete_data.set_index('datetime', inplace=True)

    # 3. Drop unnecessary columns
    complete_data.drop(columns=['time', 'DD', 'MM', 'YYYY', 'hh', 'min', 'ABSDA'], inplace=True)

    # 4. Remove non-numeric characters
    complete_data = complete_data.replace(to_replace=r'[^\d.]', value='', regex=True)

    # 5. Convert all data columns to float type
    complete_data = complete_data.astype(float)

    # 6. Replace 0.0 values with NaN
    complete_data.replace(0.0, np.nan, inplace=True)

    # Return the first few rows for verification
    complete_data.head()

    return complete_data


def read_wang(species, site, network):

    # Read ale_gage_sites.json
    with open(paths.root / "data/ale_gage_sites.json", "r") as file:
        site_info = json.load(file)

    # Read species_info.json
    with open(paths.root / "data/ale_gage_species.json", "r") as file:
        species_info = json.load(file)[species]

    site_name = site_info[site]["gcwerks_name"]
    site_code = sites_wang[site]

    data_folder = home / f"data/ale_gage_wang/{network.lower()}_data_all/complete/{site_name}"
    files = data_folder.glob(f"{site_code}-{network.lower()}*.dap")

    dfs = []

    for file in files:

        dfs.append(read_wang_file(file))

    df = pd.concat(dfs)

    # Only output one species
    df = df[[species_wang[species]]]

    # Rename species column
    df.rename(columns={species_wang[species]: "mf"}, inplace=True)

    df.index.name = "time"

    # Sort by time
    df.sort_index(inplace=True)

    # Remove duplicate indices
    df = df[~df.index.duplicated(keep='first')]

    return df


    