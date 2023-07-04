import configparser
from pathlib import Path
import xarray as xr
import json
import pandas as pd
import tarfile
import numpy as np
from datetime import datetime

from agage_archive import get_path
from agage_archive.util import is_number
from agage_archive.processing import create_dataset


class Paths():
    def __init__(self):
        """Class to store paths to data folders
        """

        # Get repository root
        self.root = get_path().parent

        # Check if config file exists
        config_file = get_path("config.ini")
        if not config_file.exists():
            raise FileNotFoundError(
                "Config file not found. Try running util.setup first")

        # Read config file
        config = configparser.ConfigParser()
        config.read(get_path("config.ini"))

        self.agage = Path(config["Paths"]["agage_path"])
        self.agage_gcmd = self.agage / "data-nc"
        self.agage_gcms = self.agage / "data-gcms-nc"
        self.ale = Path(config["Paths"]["ale_path"])
        self.gage = Path(config["Paths"]["gage_path"])

        self.output = Path(config["Paths"]["output_path"])

        # Check that data folders are there
        for pth in [self.agage_gcmd,
                    self.agage_gcms,
                    self.ale,
                    self.gage]:

            if not pth.exists():
                raise FileNotFoundError(
                    f"Folder {pth} doesn't exist")
        
        # Check that NetCDF folders contain .nc files
        for pth in [self.agage_gcmd, self.agage_gcms]:
            if not list(pth.glob("*.nc")):
                raise FileNotFoundError(
                    f"""{pth} directory doesn't
                    contain any NetCDF files"""
                )

        # Check that ALE and GAGE folder contain .C files
        for pth in [self.ale, self.gage]:
            if not list(pth.glob("*.gtar.gz")):
                raise FileNotFoundError(
                    f"""{pth} directory doesn't
                    contain any GCWerks .gtar.gz files"""
                )


paths = Paths()


def read_agage(species, site, instrument):
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

    species_search = species.lower()

    if instrument == "GCMD":
        pth = paths.agage_gcmd
    elif "GCMS" in instrument:
        pth = paths.agage_gcms

    nc_file = pth / f"AGAGE-{instrument}_{site}_{species_search}.nc"

    if not nc_file.exists():
        raise FileNotFoundError(f"Can't find file {nc_file}")

    with xr.open_dataset(nc_file) as f:
        ds = f.load()

    # Everything should have been flagged already, but just in case...
    flagged = ds.data_flag != 0
    ds.mf[flagged] = np.nan
    ds.mf_repeatability[flagged] = np.nan

    # For public files, remove flagged data and some other variables
    ds = ds.drop_vars(["data_flag",
                       "integration_flag",
                       "git_pollution_flag",
                       "met_office_baseline_flag",
                       "run_time"],
                       errors="ignore")

    return ds


def read_ale_gage(species, site, network):
    """Read GA Tech ALE/GAGE files

    Args:
        species (str): Species
        site (str): Three-letter site code
        network (str): "ALE" or "GAGE"

    Returns:
        pd.DataFrame: Pandas dataframe containing file contents
    """

    # Get data on ALE/GAGE sites
    with open(paths.root / "data/ale_gage_sites.json") as f:
        site_info = json.load(f)

    # Get species info
    with open(paths.root / "data/ale_gage_species.json") as f:
        species_info = json.load(f)[species]

    # For now, hardwire path
    #TODO: Sort out these paths
    folder = {"ALE": Path("/Users/chxmr/data/ale_gage_sio1993/ale"),
              "GAGE":  Path("/Users/chxmr/data/ale_gage_sio1993/gage")}

    pth = folder[network] / f"{site_info[site]['gcwerks_name']}_sio1993.gtar.gz"

    tar = tarfile.open(pth, "r:gz")

    dfs = []

    for member in tar.getmembers():

        # Extract tar file
        f = tar.extractfile(member)
        
        meta = f.readline().decode("ascii").strip()
        header = f.readline().decode("ascii").split()

        site_in_file = meta[:2]
        year = meta[2:4]
        month = meta[4:7]

        nspecies = len(header) - 3
        columns = header[:3]

        # Define column widths
        colspec = [3, 5, 7]
        for sp in header[3:]:
            colspec += [7, 1]
            columns += [str(sp).replace("'", ""),
                        f"{sp}_pollution"]

        # Read data
        df = pd.read_fwf(f, skiprows=0,
                        widths=colspec,
                        names=columns,
                        na_values = -99.9)

        # Some HHMM labeled as 2400, increment to following day
        midnight = df["TIME"] == 2400
        df.loc[midnight, "DA"] += 1
        df.loc[midnight, "TIME"] = 0

        # Create datetime string
        datetime = df["DA"].astype(str) + \
            f"/{month}/{year}:" + \
            df["TIME"].astype(str).str.zfill(4)

        # Convert datetime string
        # There are some strange entries in here. If we can't understand the format, reject that point.
        df.index = pd.to_datetime(datetime, format="%d/%b/%y:%H%M",
                                  errors="coerce")
        df.index = df.index.tz_localize(site_info[site]["tz"],
                                        ambiguous="NaT",
                                        nonexistent="NaT")

        dfs.append(df)

    # Concatenate monthly dataframes into single dataframe
    df_combined = pd.concat(dfs)

    # Convert to UTC
    df_combined.index = df_combined.index.tz_convert(None)

    # Sort
    df_combined.sort_index(inplace=True)

    # Drop duplicated indices
    df_combined = df_combined.loc[df_combined.index.drop_duplicates(),:]

    # Drop na indices
    df_combined = df_combined.loc[~df_combined.index.isna(),:]

    # Output one species
    df_combined = df_combined[species_info["species_name_gatech"]]

    repeatability = species_info[f"{network.lower()}_repeatability_percent"]/100.

    df_out = pd.DataFrame(index=df_combined.index,
                          data={"mf": df_combined.values.copy(),
                                "mf_repeatability": df_combined.values.copy()*repeatability})

    df_out.attrs["scale"] = species_info["scale"]
    df_out.attrs["units"] = species_info["units"]

    return create_dataset(species, site, network, df_out)


def output_dataset(ds, end_date = None):
    '''Output dataset to netCDF file

    Args:
        ds (xr.Dataset): Dataset to output
        end_date (str, optional): End date to subset to. Defaults to None.
    '''
    
    #TODO: may need to translate species
    filename = f"AGAGE-combined_{ds.attrs['site_code']}_{ds.attrs['species'].lower()}.nc"

    ds.sel(time=slice(None, end_date)).to_netcdf(paths.output / filename, mode="w", format="NETCDF4")

