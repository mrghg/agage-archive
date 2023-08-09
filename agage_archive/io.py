import xarray as xr
import json
import pandas as pd
import tarfile
import numpy as np
import pytz

from agage_archive import Paths
from agage_archive.convert import scale_convert
from agage_archive.formatting import format_species, \
    format_variables, format_attributes
from agage_archive.data_selection import read_release_schedule, data_exclude, \
    calibration_scale_default, read_data_combination
from agage_archive.definitions import instrument_type_definition


def read_agage(species, site, instrument,
               testing_path = False,
               verbose = False):
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

    paths = Paths(test=testing_path)

    species_search = format_species(species)

    gcmd_instruments = ["GCMD", "GCECD", "Picarro", "LGR"]
    gcms_instruments = ["GCMS-ADS", "GCMS-Medusa", "GCMS-MteCimone"]


    # Determine path

    pth = None

    for gcmd_instrument in gcmd_instruments:
        if gcmd_instrument in instrument:
            pth = paths.agage_gcmd
            break
    for gcms_instrument in gcms_instruments:
        if gcms_instrument in instrument:
            pth = paths.agage_gcms
            break
    
    if pth == None:
        raise ValueError(f"Instrument must be one of {gcmd_instruments} {gcms_instruments}")

    # search for netcdf files matching instrument, site and species
    nc_files = list(pth.glob(f"AGAGE-{instrument}*_{site}_{species_search}.nc"))

    if len(nc_files) == 0:
        raise FileNotFoundError(f"Can't find file AGAGE-{instrument}*_{site}_{species_search}.nc")
    elif len(nc_files) > 1:
        raise FileNotFoundError(f"Found more than one file matching AGAGE-{instrument}*_{site}_{species_search}.nc")
    else:
        nc_file = nc_files[0]

    if verbose:
        print(f"... reading {nc_file}")

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

    # Remove any excluded data
    ds = data_exclude(ds, format_species(species), site, instrument)

    # Check against release schedule
    rs = read_release_schedule(instrument,
                               species=format_species(species),
                               site=site)
    ds = ds.sel(time=slice(None, rs))


    return ds


def ale_gage_timestamp_issues(datetime, timestamp_issues):
    """Check for timestamp issues in ALE/GAGE data

    Args:
        datetime (pd.Series): Datetime series
        timestamp_issues (dict): Dictionary of timestamp issues from ale_gage_timestamp_issues.json

    Returns:
        pd.Series: Datetime series with issues fixed
    """

    if len(timestamp_issues) == 0:
        return datetime
    
    for timestamp_issue in timestamp_issues:
        if timestamp_issue in datetime.values:
            print(f"... Timestamp issue at {timestamp_issue} replacing with {timestamp_issues[timestamp_issue]}")
            datetime = datetime.replace(timestamp_issue, timestamp_issues[timestamp_issue])

    return datetime


def read_ale_gage(species, site, network,
                  testing_path = False,
                  verbose = False):
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

    paths = Paths(test=testing_path)

    # Get data on ALE/GAGE sites
    with open(paths.root / "data/ale_gage_sites.json") as f:
        site_info = json.load(f)

    # Get species info
    with open(paths.root / "data/ale_gage_species.json") as f:
        species_info = json.load(f)[format_species(species)]

    # Get Datetime issues list
    with open(paths.root / "data/ale_gage_timestamp_issues.json") as f:
        timestamp_issues = json.load(f)
        if site in timestamp_issues[network]:
            timestamp_issues = timestamp_issues[network][site]
        else:
            timestamp_issues = {}

    # Path to relevant folder
    folder = paths.__getattribute__(network.lower())

    pth = folder / f"{site_info[site]['gcwerks_name']}_sio1993.gtar.gz"

    if verbose:
        print(f"... opening {pth}")

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

            # Create datetime string. This format is a little weird, but it's an easy way to construct it
            datetime = df["DA"].astype(str).str.zfill(2) + \
                f"-{month}-{year} " + \
                df["TIME"].astype(str).str.zfill(4)
            
            # Check for datetime issues
            datetime = ale_gage_timestamp_issues(datetime, timestamp_issues)

            # Convert datetime string
            df.index = pd.to_datetime(datetime, format="%d-%b-%y %H%M")

            # Drop duplicates
            if "duplicates" in timestamp_issues:
                keep = timestamp_issues["duplicates"]
            else:
                keep = "first"
            df = df[~df.index.duplicated(keep=keep)]

            # Timestamps are local time (no daylight savings)
            tzoffset_hours = site_info[site]["tz"].split("UTC")[1]
            local_offset = pytz.FixedOffset(int(tzoffset_hours)*60)
            df.index = df.index.tz_localize(local_offset)

            dfs.append(df)

    # Concatenate monthly dataframes into single dataframe
    df_combined = pd.concat(dfs)

    # Convert to UTC
    df_combined.index = df_combined.index.tz_convert("UTC")

    # Sort
    df_combined.sort_index(inplace=True)

    # Check if there are NaN indices
    if len(df_combined.index[df_combined.index.isna()]) > 0:
        raise ValueError("NaN indices found. Check timestamp issues.")

    # Drop na indices
    df_combined = df_combined.loc[~df_combined.index.isna(),:]

    # Check if there are duplicate indices
    if len(df_combined.index) != len(df_combined.index.drop_duplicates()):
        # Find which indices are duplicated
        duplicated_indices = df_combined.index[df_combined.index.duplicated(keep=False)]
        raise ValueError(f"Duplicate indices found. Check timestamp issues: {duplicated_indices}")

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
                        species=format_species(species),
                        calibration_scale=species_info["scale"],
                        units=species_info["units"])

    ds = format_variables(ds)

    # Remove any excluded data
    ds = data_exclude(ds, format_species(species), site, network)

    # Check against release schedule
    rs = read_release_schedule(network,
                               species=format_species(species),
                               site=site)
    ds = ds.sel(time=slice(None, rs))

    return ds


def combine_datasets(species, site, 
                    scale = "default",
                    testing_path = False,
                    verbose = False):
    '''Combine ALE/GAGE/AGAGE datasets for a given species and site

    Args:
        species (str): Species
        site (str): Site
        scale (str, optional): Calibration scale. Defaults to value in scale_defaults.csv.
            If None, will attempt to leave scale unchanged.
        convert_scale (bool, optional): Convert calibration scale, or keep the same. Defaults to True.

    Returns:
        xr.Dataset: Dataset containing data
    '''

    # Read instrument dates from CSV files
    instruments = read_data_combination(format_species(species), site)

    instrument_types, instrument_number_str = instrument_type_definition()

    # Get default calibration scale, if needed
    if scale != None:
        if scale == "default":
            scale = calibration_scale_default(format_species(species))

    # Combine datasets    
    dss = []
    comments = []
    attrs = []
    instrument_rec = []
    dates_rec = []
    networks = []
    scales = []

    for instrument, date in instruments.items():

        # Read data
        if instrument in ["ALE", "GAGE"]:
            ds = read_ale_gage(species, site, instrument,
                               testing_path=testing_path,
                               verbose=verbose)
        else:
            ds = read_agage(species, site, instrument,
                            testing_path=testing_path,
                            verbose=verbose)

        # Store attributes
        attrs.append(ds.attrs)

        # Store instrument info
        instrument_rec.append({key:value for key, value in ds.attrs.items() if "instrument" in key})

        # Store comments
        comments.append(ds.attrs["comment"])

        # Store network
        networks.append(ds.attrs["network"])

        # Subset date
        ds = ds.sel(time=slice(*date))

        if len(ds.time) == 0:
            raise ValueError(f"No data retained for {species} {site} {instrument}. " + \
                             "Check dates in data_combination or omit this instrument.")
        dates_rec.append(ds.time[0].dt.strftime("%Y-%m-%d").values)

        # Convert scale
        if scale != None:
            ds = scale_convert(ds, scale)

        # Record scale
        scales.append(ds.attrs["calibration_scale"])

        # Add instrument_type to dataset as variable
        ds["instrument_type"] = xr.DataArray(np.repeat(instrument_types[instrument], len(ds.time)),
                                        dims="time", coords={"time": ds.time})

        dss.append(ds)

    # Check that we don't have different scales
    if len(set(scales)) > 1:
        error_message = "Can't combine scales that do not match. Either specify a scale, or add to scale_defaults.csv. "
        for instrument, sc in zip(instruments.keys(), scales):
            error_message += f"{instrument}:{sc}, " 
        raise ValueError(error_message)

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

    # Remove instrument_type variable if all the same    
    if len(np.unique(ds_combined.instrument_type.values)) == 1:
        ds_combined = ds_combined.drop_vars("instrument_type")

    # Update network attribute
    ds_combined.attrs["network"] = "/".join(set(networks))

    return ds_combined


def output_dataset(ds,
                   network = "AGAGE",
                   instrument = "GCMD",
                   end_date = None,
                   testing_path = False,
                   output_subpath = None,
                   verbose = False):
    '''Output dataset to netCDF file

    Args:
        ds (xr.Dataset): Dataset to output
        end_date (str, optional): End date to subset to. Defaults to None.
    '''

    paths = Paths(test=testing_path)

    if output_subpath == None:
        output_path = paths.output
    else:
        output_path = paths.output / output_subpath

    # Test if output_path exists and if not create it
    if not output_path.exists():
        output_path.mkdir(parents=True)

    # Create filename
    filename = f"{network}-{instrument}_{ds.attrs['site_code']}_{format_species(ds.attrs['species'])}.nc"

    ds_out = ds.copy(deep = True)

    # Can't have some time attributes
    if "units" in ds_out.time.attrs:
        del ds_out.time.attrs["units"]
    if "calendar" in ds_out.time.attrs:
        del ds_out.time.attrs["calendar"]

    if verbose:
        print(f"... writing {paths.output / filename}")

    # Subset time and write netCDF
    ds_out.sel(time=slice(None, end_date)).to_netcdf(output_path / filename, mode="w", format="NETCDF4")
