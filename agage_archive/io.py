import configparser
from pathlib import Path
import xarray as xr
import json
import pandas as pd

from agage_archive import get_path


class Paths():
    def __init__(self):

        # Check if config file exists
        config_file = get_path("config.ini")
        if not config_file.exists():
            raise FileNotFoundError(
                "Config file not found. Try running util.setup first")

        # Read config file
        config = configparser.ConfigParser()
        config.read(get_path("config.ini"))

        self.input_path = Path(config["Paths"]["input_path"])
        self.output_path = Path(config["Paths"]["output_path"])

        # Data folders expected within input path
        self.data_suffix = "data-nc"
        self.data_gcms_suffix = "data-gcms-nc"
        self.data_ale_suffix = "data-ale"
        self.data_gage_suffix = "data-gage"

        # Check that data folders are there
        for suffix in [self.data_suffix,
                       self.data_gcms_suffix,
                       self.data_ale_suffix,
                       self.data_gage_suffix]:

            if not (self.input_path / suffix).exists():
                raise FileNotFoundError(
                    f"Data folder must contain {suffix} folder")
        
        # Check that NetCDF folders contain .nc files
        for suffix in [self.data_suffix, self.data_gcms_suffix]:
            if not list((self.input_path / suffix).glob("*.nc")):
                raise FileNotFoundError(
                    f"""{suffix} directory doesn't
                    contain any NetCDF files"""
                )

        # Check that ALE and GAGE folder contain .C files
        for suffix in [self.data_ale_suffix, self.data_gage_suffix]:
            if not list((self.input_path / suffix).glob("*.C")):
                raise FileNotFoundError(
                    f"""{suffix} directory doesn't
                    contain any GCWerks .C files"""
                )

paths = Paths()


def read_nc(species, site, instrument):

    nc_file = paths.input_path / paths.data_suffix
    nc_file = nc_file / f"AGAGE-{instrument}_{site}_{species}.nc"

    if not nc_file.exists():
        raise FileNotFoundError(f"Can't find file {nc_file}")

    with xr.open_dataset(nc_file) as f:
        ds = f.load()

    return ds


def read_c(site, network):
    """Read .C files containing ALE/GAGE data

    Args:
        site (str): site code

    Returns:
        pd.Dataframe: Pandas Dataframe containing concatenated contents of .C files
    """

    # Get data on ALE/GAGE sites
    with open(get_path().parent / "data/ale_gage_sites.json") as f:
        site_info = json.load(f)

    # Find, open and concatenate files for site
    c_files = []
    search_string = f"{site_info[site]['gcwerks_name']}_{network}.*.C"
    suffix = getattr(paths, f"data_{network}_suffix")
    c_files += (paths.input_path / suffix).glob(search_string)

    dfs = []
    for c_file in c_files:
        dfs.append(pd.read_fwf(c_file,
            skiprows=4,
            colspecs="infer",
            parse_dates={"datetime": ["yyyy", "mm", "dd", "hh", "mi"]}))
    df = pd.concat(dfs)

    # Create time index
    df.index = pd.to_datetime(df["datetime"], format="%Y %m %d %H %M")
    df.index.name = None

    # Rename and drop columns
    columns = []
    for ci, col in enumerate(df.columns):
        if "Flag" in col:
            columns.append(f"{df.columns[ci-1]}_flag")
        else:
            columns.append(col)
    df.columns = columns

    df.drop(columns=["Year", "Inlet", "Standard", "datetime"], inplace=True)

    return df.sort_index()
