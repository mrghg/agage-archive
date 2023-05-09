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
    """Read GCWerks netCDF files

    Args:
        species (str): Species
        site (str): Site code
        instrument (str): Instrument

    Raises:
        FileNotFoundError: Can't find netCDF file

    Returns:
        xarray.Dataset: Contents of netCDF file
    """

    nc_file = paths.input_path / paths.data_suffix
    nc_file = nc_file / f"AGAGE-{instrument}_{site}_{species}.nc"

    if not nc_file.exists():
        raise FileNotFoundError(f"Can't find file {nc_file}")

    with xr.open_dataset(nc_file) as f:
        ds = f.load()

    return ds


def read_c(species, site, network):
    """Read .C files containing ALE/GAGE data

    Args:
        site (str): site code
        network (str): network, either ALE or GAGE

    Returns:
        pd.Dataframe: Pandas Dataframe containing concatenated contents of .C files
    """

    # Get data on ALE/GAGE sites
    with open(get_path().parent / "data/ale_gage_sites.json") as f:
        site_info = json.load(f)

    # Find, open and concatenate files for site
    c_files = []
    search_string = f"{site_info[site]['gcwerks_name']}_{network.lower()}.*.C"
    suffix = getattr(paths, f"data_{network.lower()}_suffix")
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

    df = df.rename(columns={species: "mf"})
    #TODO: This needs to be done properly!
    df[f"mf_uncertainty"] = df["mf"]*0.02

    return df[["mf", "mf_uncertainty"]].sort_index()
#    return df[[species, f"{species}_uncertainty"]].sort_index()


def combine_datasets(species, site):

    # Get instructions on how to combine datasets
    with open(get_path().parent / "data/data_selector.json") as f:
        data_selector = json.load(f)

    instruments = [["1970-01-01", "Medusa"]]

    if site in data_selector:
        if species in data_selector[site]:
            instruments = data_selector[site][species]
    
    dfs = []
    for i, (date, instrument) in enumerate(instruments):
        if instrument in ["ALE", "GAGE"]:
            df = read_c(species, site, instrument)
        else:
            df = read_nc(species, site, instrument)

        # Apply start date
        df = df.loc[date:, :]

        # Cut off previous instrument timeseries
        if i > 0:
            dfs[i-1] = dfs[i-1].loc[:date, :]

        dfs.append(df)

    return dfs


