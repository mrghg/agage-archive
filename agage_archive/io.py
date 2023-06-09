import configparser
from pathlib import Path
import xarray as xr
import json
import pandas as pd
import tarfile
import numpy as np

from agage_archive import Paths
from agage_archive.processing import create_dataset, global_attributes_instrument, \
    global_attributes_combine_instruments, scale_convert, format_dataset,\
    format_species

instrument_number = {"ALE": 0,
                    "GAGE": 1,
                    "GCMD": 2,
                    "GCMS-ADS": 3,
                    "GCMS-Medusa": 4}

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

    # Read sampling time
    if "sampling_time_seconds" in ds.time.attrs:
        sampling_time = int(ds.time.attrs["sampling_time_seconds"])
    else:
        # GCMD files don't have sampling time in the file
        # assume it's 30s (Peter Salameh, pers. comm., 2023-07-06)
        sampling_time = 30

    # Add sampling time to variables
    ds["sampling_time"] = xr.DataArray(np.ones(len(ds.time)).astype(np.int16)*sampling_time,
                                        coords={"time": ds.time},
                                        attrs={"units": "s",
                                            "long_name": "sampling_time",
                                            "comment": "Sampling period in seconds"})
    ds["sampling_time"].encoding = {"dtype": "int16"}

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

    # Add instrument global attributes, if needed
    ds = global_attributes_instrument(ds, instrument)

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

    return create_dataset(df_out, species, site, network, f"{network} GCMD")


def combine_datasets(species, site, scale = "SIO-05"):
    '''Combine ALE/GAGE/AGAGE datasets for a given species and site

    Args:
        species (str): Species
        site (str): Site
        scale (str, optional): Calibration scale. Defaults to "SIO-05".
            If None, no scale conversion is attempted

    Returns:
        xr.Dataset: Dataset containing data
    '''

    # Get instructions on how to combine datasets
    with open(paths.root / "data/data_selector.json") as f:
        data_selector = json.load(f)
    
    # Set default to Medusa
    instruments = {"GCMS-Medusa": ["", ""]}

    # Read instruments from JSON file
    if site in data_selector:
        if species in data_selector[site]:
            instruments = data_selector[site][species]
    
    dss = []
    comments = []
    attrs = []
    instrument_rec = []
    dates_rec = []

    for instrument, date in instruments.items():

        if instrument in ["ALE", "GAGE"]:
            ds = read_ale_gage(species, site, instrument)
        else:
            ds = read_agage(species, site, instrument)

        attrs.append(ds.attrs)

        instrument_rec.append({key:value for key, value in ds.attrs.items() if "instrument" in key})

        comments.append(ds.attrs["comment"])

        # Subset date
        date = [None if d == "" else d for d in date]
        ds = ds.sel(time=slice(*date))

        dates_rec.append(ds.time[0].dt.strftime("%Y-%m-%d").values)

        # Convert scale
        if scale != None:
            ds = scale_convert(ds, scale)

        # Add instrument_type to dataset as variable
        ds["instrument_type"] = xr.DataArray(np.repeat(instrument_number[instrument], len(ds.time)),
                                        dims="time", coords={"time": ds.time})
        ds.instrument_type.encoding = {"dtype": "int8"}
        ds.instrument_type.attrs["long_name"] = "ALE/GAGE/AGAGE instrument type"
        ds.instrument_type.attrs["comment"] = "0 = GC multi-detector (GCMD) from the ALE project; " + \
                                        "1 = GCMD from the GAGE project; " + \
                                        "2, 3, 4 are GCMD, GCMS-ADS or GCMS-Medusa instruments from AGAGE, respectively"
        ds.instrument_type.attrs["units"] = ""

        dss.append(ds)

    ds_combined = xr.concat(dss, dim="time")

    # Sort by time
    ds_combined = ds_combined.sortby("time")

    # Add details on instruments to global attributes
    ds_combined = global_attributes_combine_instruments(ds_combined,
                                                        instrument_rec)

    # Extend comment attribute describing all datasets
    if len(comments) > 1:
        ds_combined.attrs["comment"] = "Combined AGAGE/GAGE/ALE dataset from the following individual sources: ---- " + \
            "|---| ".join(comments)
    else:
        ds_combined.attrs["comment"] = comments[0]

    # Add site code
    ds_combined.attrs["site_code"] = site.upper()

    return ds_combined


def output_dataset(ds,
                   network = "AGAGE",
                   instrument = "GCMD",
                   end_date = None):
    '''Output dataset to netCDF file

    Args:
        ds (xr.Dataset): Dataset to output
        end_date (str, optional): End date to subset to. Defaults to None.
    '''

    #TODO: is this the best place to call this?    
    ds = format_dataset(ds)

    #TODO: maybe add network and instrument to attributes so it doesn't have to be an input?
    # NOTE: having said that instrument may not be reliable, as it's used inconsistently in AGAGE files
    filename = f"{network}-{instrument}_{ds.attrs['site_code']}_{format_species(ds.attrs['species'])}.nc"

    ds.sel(time=slice(None, end_date)).to_netcdf(paths.output / filename, mode="w", format="NETCDF4")

