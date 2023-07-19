import configparser
from pathlib import Path
import xarray as xr
import json
import pandas as pd
import tarfile
import numpy as np

from agage_archive import Paths
from agage_archive.processing import format_species, \
    format_variables, format_attributes, scale_convert


# Define instrument numbers
# NOTE: If this changes, make sure you update the instrument_type comment in variables.json
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
        sampling_period = int(ds.time.attrs["sampling_time_seconds"])
    else:
        # GCMD files don't have sampling time in the file
        # assume it's 30s (Peter Salameh, pers. comm., 2023-07-06)
        sampling_period = 30

    # Add sampling time to variables
    ds["sampling_period"] = xr.DataArray(np.ones(len(ds.time)).astype(np.int16)*sampling_period,
                                        coords={"time": ds.time})

    # Everything should have been flagged already, but just in case...
    flagged = ds.data_flag != 0
    ds.mf[flagged] = np.nan
    ds.mf_repeatability[flagged] = np.nan

    ds.attrs["site_code"] = site.upper()

    ds = format_variables(ds)

    # If no instrument attributes are present, add them using format_attributes
    if "instrument" not in ds.attrs:
        instruments = [{"instrument": instrument}]
    else:
        instruments = []

    ds = format_attributes(ds,
                        instruments=instruments,
                        species=species)

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

    if network not in ["ALE", "GAGE"]:
        raise ValueError("network must be ALE or GAGE")

    # Get data on ALE/GAGE sites
    with open(paths.root / "data/ale_gage_sites.json") as f:
        site_info = json.load(f)

    # Get species info
    with open(paths.root / "data/ale_gage_species.json") as f:
        species_info = json.load(f)[species]

    # For now, hardwire path
    folder = paths.__getattribute__(network.lower())

    pth = folder / f"{site_info[site]['gcwerks_name']}_sio1993.gtar.gz"

    with tarfile.open(pth, "r:gz") as tar:

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

    # Estimate of repeatability
    repeatability = species_info[f"{network.lower()}_repeatability_percent"]/100.

    nt = len(df_combined.index)

    # Create xarray dataset
    ds = xr.Dataset(data_vars={"mf": ("time", df_combined.values.copy()),
                            "mf_repeatability": ("time", df_combined.values.copy()*repeatability),
                            "inlet_height": ("time", np.repeat(site_info[site]["inlet_height"], nt)),
                            "sampling_period": ("time", np.repeat(30, nt)),
                            },
                    coords={"time": df_combined.index.copy()})

    # Global attributes
    comment = f"{network} {species} data from {site_info[site]['station_long_name']}. " + \
        "This data was originally processed by Derek Cunnold, Georgia Institute of Technology, " + \
        "from the original files and has now been reprocessed into netCDF format."

    ds.attrs = {"comment": comment,
                "data_owner_email": site_info[site]["data_owner_email"],
                "data_owner": site_info[site]["data_owner"],
                "station_long_name": site_info[site]["station_long_name"],
                "inlet_base_elevation_masl": site_info[site]["inlet_base_elevation_masl"],
                "inlet_latitude": site_info[site]["latitude"],
                "inlet_longitude": site_info[site]["longitude"],
                "inlet_comment": "",
                "network": network,
                "site_code": site}

    ds = format_attributes(ds,
                        instruments=[{"instrument": f"{network.upper()}_GCMD"}],
                        species=species,
                        calibration_scale=species_info["scale"],
                        units=species_info["units"])

    ds = format_variables(ds)

    return ds


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

        dss.append(ds)

    ds_combined = xr.concat(dss, dim="time")

    # Sort by time
    ds_combined = ds_combined.sortby("time")

    # Add details on instruments to global attributes
    ds_combined = format_attributes(ds_combined,
                                instrument_rec)

    # Extend comment attribute describing all datasets
    if len(comments) > 1:
        comment_str = "Combined AGAGE/GAGE/ALE dataset from the following individual sources:\n"
        for i, comment in enumerate(comments):
            comment_str += f"{i}) {comment}\n"
    else:
        comment_str = comments[0]

    ds_combined.attrs["comment"] = comment_str

    # Format variables
    ds_combined = format_variables(ds_combined)

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

    filename = f"{network}-{instrument}_{ds.attrs['site_code']}_{format_species(ds.attrs['species'])}.nc"

    ds_out = ds.copy(deep = True)

    if "units" in ds_out.time.attrs:
        del ds_out.time.attrs["units"]
        del ds_out.time.attrs["calendar"]


    # Remove instrument_type variable if all the same    
    if len(np.unique(ds_out.instrument_type.values)) == 1:
        ds_out = ds_out.drop_vars("instrument_type")

    ds_out.sel(time=slice(None, end_date)).to_netcdf(paths.output / filename, mode="w", format="NETCDF4")
