import xarray as xr
import json
import pandas as pd
import numpy as np

from agage_archive import Paths, open_data_file, data_file_list
from agage_archive.convert import scale_convert
from agage_archive.formatting import format_species, \
    format_variables, format_attributes
from agage_archive.data_selection import read_release_schedule, read_data_exclude, \
    read_data_combination
from agage_archive.definitions import instrument_type_definition
from agage_archive.util import tz_local_to_utc


gcwerks_species = {"c2f6": "pfc-116",
                   "c3f8": "pfc-218",
                   "c4f8": "pfc-318"}
    

def read_nc(network, species, site, instrument,
            verbose = False,
            data_exclude = True,
            scale = "default"):
    """Read GCWerks netCDF files

    Args:
        network (str): Network, e.g., "agage"
        species (str): Species
        site (str): Site code
        instrument (str): Instrument
        verbose (bool, optional): Print verbose output. Defaults to False.
        data_exclude (bool, optional): Exclude data based on data_exclude.xlsx. Defaults to True.
        scale (str, optional): Scale to convert to. Defaults to "default". If None, will keep original scale.

    Raises:
        FileNotFoundError: Can't find netCDF file

    Returns:
        xarray.Dataset: Contents of netCDF file
    """

    paths = Paths(network)

    species_search = format_species(species)
    if species_search in gcwerks_species:
        species_search = gcwerks_species[species_search]

    gcmd_instruments = ["GCMD", "GCECD", "Picarro", "LGR"]
    gcms_instruments = ["GCMS-ADS", "GCMS-Medusa", "GCMS-MteCimone"]

    # Determine sub-path within data directory
    sub_path = None

    for gcmd_instrument in gcmd_instruments:
        if gcmd_instrument in instrument:
            sub_path = paths.agage_md_path
            break
    for gcms_instrument in gcms_instruments:
        if gcms_instrument in instrument:
            sub_path = paths.agage_gcms_path
            break
    
    if sub_path == None:
        raise ValueError(f"Instrument must be one of {gcmd_instruments} {gcms_instruments}")

    # search for netcdf files matching instrument, site and species
    nc_files = data_file_list(network, sub_path, f"*-{instrument}*_{site}_{species_search}.nc")[2]

    if len(nc_files) == 0:
        raise FileNotFoundError(f"Can't find file matching *-{instrument}*_{site}_{species_search}.nc in data/{network}/{sub_path}")
    elif len(nc_files) > 1:
        raise FileNotFoundError(f"Found more than one file matching *-{instrument}*_{site}_{species_search}.nc in data/{network}/{sub_path}")
    else:
        nc_file = nc_files[0]

    if verbose:
        print(f"... reading {nc_file}")

    # Read netCDF file
    with open_data_file(nc_file, network, sub_path=sub_path, verbose=True) as f:
        with xr.open_dataset(f, engine="h5netcdf") as ds_file:
            ds = ds_file.load()

    # Read sampling time
    if "sampling_time_seconds" in ds.time.attrs:
        sampling_period = int(ds.time.attrs["sampling_time_seconds"])
    else:
        # GCMD files don't have sampling time in the file
        # assume it's 1s (Peter Salameh, pers. comm., 2023-07-06)
        sampling_period = 1

    # Add sampling time to variables
    ds["sampling_period"] = xr.DataArray(np.ones(len(ds.time)).astype(np.int16)*sampling_period,
                                        coords={"time": ds.time})

    # Everything should have been flagged in the AGAGE files already, but just in case...
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
                        network=network,
                        species=species)

    # Remove any excluded data
    if data_exclude:
        ds = read_data_exclude(ds, format_species(species), site, instrument)

    # Check against release schedule
    rs = read_release_schedule(network, 
                               instrument,
                               species=format_species(species),
                               site=site)
    ds = ds.sel(time=slice(None, rs))

    # Convert scale, if needed
    ds = scale_convert(ds, scale)

    return ds


def ale_gage_timestamp_issues(datetime, timestamp_issues,
                              verbose = True):
    """Check for timestamp issues in ALE/GAGE data

    Args:
        datetime (pd.Series): Datetime series
        timestamp_issues (dict): Dictionary of timestamp issues from ale_gage_timestamp_issues.json
        verbose (bool, optional): Print verbose output. Defaults to False.
        
    Returns:
        pd.Series: Datetime series with issues fixed
    """

    if len(timestamp_issues) == 0:
        return datetime
    
    for timestamp_issue in timestamp_issues:
        if timestamp_issue in datetime.values:
            if verbose:
                print(f"... Timestamp issue at {timestamp_issue} replacing with {timestamp_issues[timestamp_issue]}")
            datetime = datetime.replace(timestamp_issue, timestamp_issues[timestamp_issue])

    return datetime


def read_ale_gage(network, species, site, instrument,
                  verbose = True,
                  utc = True,
                  data_exclude = True,
                  scale = "default"):
    """Read GA Tech ALE/GAGE files, process and clean

    Args:
        network (str): Network. Can only be "agage" or "agage_test"
        species (str): Species
        site (str): Three-letter site code
        instrument (str): "ALE" or "GAGE"
        verbose (bool, optional): Print verbose output. Defaults to False.
        utc (bool, optional): Convert to UTC. Defaults to True.
        data_exclude (bool, optional): Exclude data based on data_exclude.xlsx. Defaults to True. 
            utc must also be true, as timestamps are in UTC.
        scale (str, optional): Calibration scale. Defaults to None, which means no conversion is attempted.
            Set to "default" to use value in scale_defaults.csv.

    Returns:
        pd.DataFrame: Pandas dataframe containing file contents
    """

    if network not in ["agage", "agage_test"]:
        raise ValueError("network must be agage or agage_test")
    
    if instrument not in ["ALE", "GAGE"]:
        raise ValueError("instrument must be ALE or GAGE")

    paths = Paths(network)

    # Get data on ALE/GAGE sites
    with open_data_file("ale_gage_sites.json", network = network) as f:
        site_info = json.load(f)

    # Get species info
    with open_data_file("ale_gage_species.json", network = network) as f:
        species_info = json.load(f)[format_species(species)]

    # Get Datetime issues list
    with open_data_file("ale_gage_timestamp_issues.json", network = network) as f:
        timestamp_issues = json.load(f)
        if site in timestamp_issues[instrument]:
            timestamp_issues = timestamp_issues[instrument][site]
        else:
            timestamp_issues = {}

    # Path to relevant sub-folder
    folder = paths.__getattribute__(f"{instrument.lower()}_path")

    with open_data_file(f"{site_info[site]['gcwerks_name']}_sio1993.gtar.gz",
                        network = network,
                        sub_path = folder,
                        verbose=True) as tar:

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
            datetime = ale_gage_timestamp_issues(datetime, timestamp_issues,
                                                 verbose=verbose)

            # Convert datetime string
            df.index = pd.to_datetime(datetime, format="%d-%b-%y %H%M")

            # Drop duplicates
            if "duplicates" in timestamp_issues:
                keep = timestamp_issues["duplicates"]
            else:
                keep = "first"
            df = df[~df.index.duplicated(keep=keep)]

            # Timestamps are local time (no daylight savings)
            if utc:
                df.index = tz_local_to_utc(df.index, network, site)

            dfs.append(df)

    # Concatenate monthly dataframes into single dataframe
    df_combined = pd.concat(dfs)

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
    repeatability = species_info[f"{instrument.lower()}_repeatability_percent"]/100.

    nt = len(df_combined.index)

    # Create xarray dataset
    ds = xr.Dataset(data_vars={"mf": ("time", df_combined.values.copy()),
                            "mf_repeatability": ("time", df_combined.values.copy()*repeatability),
                            "inlet_height": ("time", np.repeat(site_info[site]["inlet_height"], nt)),
                            "sampling_period": ("time", np.repeat(30, nt)),
                            },
                    coords={"time": df_combined.index.copy()})

    # Global attributes
    comment = f"{instrument} {species} data from {site_info[site]['station_long_name']}. " + \
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
                "site_code": site}

    ds = format_attributes(ds,
                        instruments=[{"instrument": f"{instrument.upper()}_GCMD"}],
                        network=network,
                        species=format_species(species),
                        calibration_scale=species_info["scale"],
                        units=species_info["units"])

    ds = format_variables(ds)

    # Remove any excluded data. Only do this if time is UTC, otherwise it won't be in the file
    if data_exclude:
        if not utc:
            raise ValueError("Can't exclude data if time is not UTC")
        ds = read_data_exclude(ds, format_species(species), site, instrument)

    # Check against release schedule
    rs = read_release_schedule(network, instrument,
                               species=format_species(species),
                               site=site)
    ds = ds.sel(time=slice(None, rs))

    # Convert scale, if needed
    ds = scale_convert(ds, scale)

    return ds


def combine_datasets(network, species, site, 
                    scale = "default",
                    verbose = True):
    '''Combine ALE/GAGE/AGAGE datasets for a given species and site

    Args:
        network (str): Network
        species (str): Species
        site (str): Site
        scale (str, optional): Calibration scale. Defaults to value in scale_defaults.csv.
            If None, will attempt to leave scale unchanged.
        verbose (bool, optional): Print verbose output. Defaults to False.

    Returns:
        xr.Dataset: Dataset containing data
    '''

    # Read instrument dates from CSV files
    instruments = read_data_combination(network, format_species(species), site)

    instrument_types, instrument_number_str = instrument_type_definition()

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
            ds = read_ale_gage(network, species, site, instrument,
                               verbose=verbose,
                               scale=scale)
        else:
            ds = read_nc(network, species, site, instrument,
                        verbose=verbose,
                        scale=scale)

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


def output_dataset(network, ds,
                   instrument = "GCMD",
                   end_date = None,
                   verbose = False):
    '''Output dataset to netCDF file

    Args:
        ds (xr.Dataset): Dataset to output
        network (str): Network
        instrument (str, optional): Instrument. Defaults to "GCMD".
        end_date (str, optional): End date to subset to. Defaults to None.
    '''

    paths = Paths(network)

    output_path = paths.data / network / paths.output_path

    # Test if output_path exists and if not create it
    if not output_path.exists():
        output_path.mkdir(parents=True)

    # Create filename
    filename = f"{network.upper()}-{instrument}_{ds.attrs['site_code']}_{format_species(ds.attrs['species'])}.nc"

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
