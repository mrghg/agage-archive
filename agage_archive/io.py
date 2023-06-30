import configparser
from pathlib import Path
import xarray as xr
import json
import pandas as pd
import tarfile
import numpy as np
from datetime import datetime


from agage_archive import get_path


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


def scale_convert(ds, scale_new):
    """Convert mole fraction from one scale to another

    Args:
        species (str): Species
        scale_original (str): Original scale
        scale_new (str): Output scale
        t (pd.Timestamp): Timestamp
        mf (float): Mole fraction
    
    Returns:
        ndarray,float: Mole fraction in new scale
    """

    def n2o_scale_function(time):
        """Function to apply to N2O mole fractions to convert from SIO-93 to SIO-98

        Args:
            time (pd.Timestamp): Timestamp
            mf (ndarray): Mole fractions

        Returns:
            ndarray: Mole fractions adjusted for time-variation between scales (excluding factor)
        """

        # calculate days elapsed since 19th August 1977
        days_since_ale_start = (time - \
                                pd.Timestamp("1978-03-02")).dt.days.values

        a0=1.00664
        a1=-0.56994e-3
        a2=-0.65398e-3
        a3=0.13083e-3
        a4=-0.20742e-4

        t = (days_since_ale_start-3227.)/365.25
        f = 1./(a0 + a1*t + a2*t**2 + a3*t**3 + a4*t**4)

        # Apply f to mf only between 1st May 1984 and 31st March 1990
        f_out = np.ones_like(f).astype(float)
        idx = (days_since_ale_start>=2252) & (days_since_ale_start<=4412)
        f_out[idx] = f[idx]
        #mf[idx] = mf[idx] * f[idx]

        return f_out
    
    # Check if scales are the same
    if ds.attrs["calibration_scale"] == scale_new:
        return ds
    else:
        scale_original = ds.attrs["calibration_scale"]
    
    species = ds.attrs["species"]

    # Make a deep copy of the dataset
    ds_out = ds.copy(deep=True)

    # Read scale conversion factors
    scale_converter = pd.read_csv(paths.root / "data/scale_convert.csv",
                                  index_col="Species")

    scale_numerator = [ratio.split("/")[0] for ratio in scale_converter.columns]
    scale_denominator = [ratio.split("/")[1] for ratio in scale_converter.columns]

    # Check for duplicates in numerator or denominator (can't handle this yet)
    if (len(set(scale_denominator)) != len(scale_denominator)) or \
        (len(set(scale_numerator)) != len(scale_numerator)):
        raise NotImplementedError("Can't deal with multiple factors for same scale at the moment")

    # Find chain of ratios to apply (start from end and work backwards)
    columns = [scale_numerator.index(scale_new)]
    while scale_denominator[columns[-1]] != scale_original:
        columns.append(scale_numerator.index(scale_denominator[columns[-1]]))

    # Now reverse to propagate forwards
    columns = columns[::-1]    

    # Apply scale conversion factors
    for column in columns:
        column_name = scale_converter.columns[column]

        if species.lower() == "n2o" and column_name == "SIO-98/SIO-93":
            # Apply time-varying factor to N2O (scale factor is not included)
            ds_out.mf.values *= n2o_scale_function(ds.time.to_series())

        ds_out.mf.values *= scale_converter.loc[species, column_name]

    # Update attributes
    ds_out.attrs["calibration_scale"] = scale_new

    return ds_out


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


def create_dataset(species, site, network, df):
    '''Create xarray dataset from pandas dataframe
    Need the following attributes:
        "comment"
        "data_owner_email"
        "station_long_name"
        "inlet_base_elevation"
        "inlet_latitude"
        "inlet_longitude"
        "inlet_comment"
        "data_dir"
        "species"
        "calibration_scale"
        "units"
        "file_created"
    and the following variables:
        "inlet_height"
        "mf"
        "mf_repeatability"
        "data_flag"
        "integration_flag"

    Args:
        species (str): Species
        site (str): Site
        network (str): Network
        df (pd.DataFrame): Dataframe containing data

    Returns:
        xr.Dataset: Dataset containing data
    '''

    units = df.attrs["units"]
    scale = df.attrs["scale"]
    github_url = "https://github.com/mrghg/agage-archive"

    # Get data on ALE/GAGE sites
    with open(paths.root / "data/ale_gage_sites.json") as f:
        site_info = json.load(f)

    # Get species info
    with open(paths.root / "data/ale_gage_species.json") as f:
        species_info = json.load(f)[species]

    nt = len(df.index)
    inlet_height = np.repeat(site_info[site]["inlet_height"], nt)
    data_flag = np.repeat(0, nt)
    integration_flag = np.repeat(0, nt)

    # Create xarray dataset
    ds = xr.Dataset(data_vars={"mf": ("time", df["mf"].values.copy()),
                            "mf_repeatability": ("time", df["mf_repeatability"].values.copy()),
                            "inlet_height": ("time", inlet_height),
                            },
                    coords={"time": df.index})
    
    # Variable encoding
    ds.mf.encoding = {"dtype": "float32"}
    ds.mf_repeatability.encoding = {"dtype": "float32"}
    ds.inlet_height.encoding = {"dtype": "int8"}
    ds.time.encoding = {"units": f"seconds since 1970-01-01 00:00:00"}

    # Variable attributes
    attributes = {"time": {"label": "centre",
                        "standard_name": "time",
                        "comment": "Time stamp corresponds to middle of sampling period. " + \
                        "Time since midnight UTC of reference date."},
                    "mf": {"units": units,
                        "scale": scale,
                        "long_name": f"mole_fraction_of_{species.lower()}_in_air",
                        "comment": f"Mole fraction of {species.lower()} in dry air"},
                    "mf_repeatability": {"units": units,
                        "long_name": f"mole_fraction_of_{species.lower()}_in_air_repeatability",
                        "comment": f"Repeatability of {species.lower()} measurements in dry air"},
                    "inlet_height": {"units": "m",
                        "long_name": f"inlet_height",
                        "comment": f"Height of inlet above inlet_base_elevation_masl"},
                    }

    for var in attributes:
        ds[var].attrs.update(attributes[var])

    # Global attributes
    comment = f"{network} {species} data from {site_info[site]['station_long_name']}. " + \
        "This data was processed by Derek Cunnold, Georgia Institute of Technology, " + \
        "from the original files and has now been reprocessed into netCDF format " + \
        f"using code at {github_url}."

    global_attributes = {"comment": comment,
                        "data_owner_email": site_info[site]["data_owner_email"],
                        "data_owner": site_info[site]["data_owner"],
                        "station_long_name": site_info[site]["station_long_name"],
                        "inlet_base_elevation_masl": site_info[site]["inlet_base_elevation_masl"],
                        "inlet_latitude": site_info[site]["latitude"],
                        "inlet_longitude": site_info[site]["longitude"],
                        "inlet_comment": "",
                        "data_dir": "",
                        "species": species.lower(),
                        "calibration_scale": scale,
                        "units": units,
                        "file_created": datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

    ds.attrs.update(global_attributes)    

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

    instrument_number = {"ALE": 0,
                         "GAGE": 1,
                         "GCMD": 2,
                         "GCMS-ADS": 3,
                         "GCMS-Medusa": 4}
    
    # Set default to Medusa
    instruments = {"GCMS-Medusa": ["", ""]}

    # Read instruments from JSON file
    if site in data_selector:
        if species in data_selector[site]:
            instruments = data_selector[site][species]
    
    dss = []
    comments = []

    for instrument, date in instruments.items():

        if instrument in ["ALE", "GAGE"]:
            ds = read_ale_gage(species, site, instrument)
        else:
            ds = read_agage(species, site, instrument)

        comments.append(ds.attrs["comment"])

        # Subset date
        date = [None if d == "" else d for d in date]
        ds = ds.sel(time=slice(*date))

        # Convert scale
        if scale != None:
            ds = scale_convert(ds, scale)

        # Add instrument to dataset as variable
        ds["instrument"] = xr.DataArray(np.repeat(instrument_number[instrument], len(ds.time)),
                                        dims="time", coords={"time": ds.time})
        ds.instrument.encoding = {"dtype": "int8"}
        ds.instrument.attrs["long_name"] = "ALE/GAGE/AGAGE instrument"
        ds.instrument.attrs["comment"] = "0 = GC multi-detector (GCMD) from the ALE project; " + \
                                        "1 = GCMD from the GAGE project; " + \
                                        "2, 3, 4 are GCMD, GCMS-ADS or GCMS-Medusa instruments from AGAGE, respectively"
        ds.instrument.attrs["units"] = ""

        dss.append(ds)

    ds_combined = xr.concat(dss, dim="time")

    # Sort by time
    ds_combined = ds_combined.sortby("time")

    # Extend comment attribute describing all datasets
    ds_combined.attrs["comment"] = "Combined AGAGE/GAGE/ALE dataset combined from the following individual sources: ---- " + \
        "|---| ".join(comments)
    
    # Add site code
    ds_combined.attrs["site_code"] = site.upper()

    return ds_combined


def output_dataset(ds, end_date = None):
    '''Output dataset to netCDF file

    Args:
        ds (xr.Dataset): Dataset to output
        end_date (str, optional): End date to subset to. Defaults to None.
    '''
    
    #TODO: may need to translate species
    filename = f"AGAGE-combined_{ds.attrs['site_code']}_{ds.attrs['species'].lower()}.nc"

    ds.sel(time=slice(None, end_date)).to_netcdf(paths.output / filename, mode="w", format="NETCDF4")

