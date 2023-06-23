import configparser
from pathlib import Path
import xarray as xr
import json
import pandas as pd
import tarfile
from tzwhere import tzwhere
import numpy as np

from agage_archive import get_path


class Paths():
    def __init__(self):

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

        self.input_path = Path(config["Paths"]["input_path"])
        self.output_path = Path(config["Paths"]["output_path"])

        # Data folders expected within input path
        self.data_suffix = "data-nc"
        self.data_gcms_suffix = "data-gcms-nc"
        self.data_ale_suffix = "data-ale"
        self.data_gage_suffix = "data-gage"

        # # Check that data folders are there
        # for suffix in [self.data_suffix,
        #                self.data_gcms_suffix,
        #                self.data_ale_suffix,
        #                self.data_gage_suffix]:

        #     if not (self.input_path / suffix).exists():
        #         raise FileNotFoundError(
        #             f"Data folder must contain {suffix} folder")
        
        # Check that NetCDF folders contain .nc files
        for suffix in [self.data_suffix, self.data_gcms_suffix]:
            if not list((self.input_path / suffix).glob("*.nc")):
                raise FileNotFoundError(
                    f"""{suffix} directory doesn't
                    contain any NetCDF files"""
                )

        # # Check that ALE and GAGE folder contain .C files
        # for suffix in [self.data_ale_suffix, self.data_gage_suffix]:
        #     if not list((self.input_path / suffix).glob("*.C")):
        #         raise FileNotFoundError(
        #             f"""{suffix} directory doesn't
        #             contain any GCWerks .C files"""
        #         )


paths = Paths()


def scale_convert(species, scale_original, scale_new, t, mf):
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

    def n2o_scale_function(time, mf):
        """Function to apply to N2O mole fractions to convert from SIO-93 to SIO-98

        Args:
            time (pd.Timestamp): Timestamp
            mf (ndarray): Mole fractions

        Returns:
            ndarray: Mole fractions adjusted for time-variation between scales (excluding factor)
        """

        # ensure mf is array
        mf = np.atleast_1d(mf)

        # calculate days elapsed since 19th August 1977
        days_since_ale_start = np.atleast_1d((time - pd.Timestamp("1978-03-02")).days)

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
    if scale_original == scale_new:
        return mf

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
            mf *= n2o_scale_function(t, mf)

        mf *= scale_converter.loc[species, column_name]

    return mf


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


def read_ale_gage(species, site, network):
    """Read GA Tech ALE/GAGE files

    Args:
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
        for species in header[3:]:
            colspec += [7, 1]
            columns += [str(species).replace("'", ""),
                        f"{species}_pollution"]

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
        df.index = pd.to_datetime(datetime, format="%d/%b/%y:%H%M", errors="coerce")
        df.index = df.index.tz_localize(site_info[site]["tz"])

        dfs.append(df)

    # Concatenate monthly dataframes into single dataframe
    df_combined = pd.concat(dfs)

    # Convert to UTC
    df_combined.index = df_combined.index.tz_convert(None)

    # Sort
    df_combined.sort_index(inplace=True)

    # Output one species
    df_combined = df_combined[species_info["species_name_gatech"]]

    repeatability = species_info[f"{network.lower()}_repeatability_percent"]/100.

    df_out = pd.DataFrame(index=df_combined.index,
                          data={"mf": df_combined.values,
                                "mf_repeatability": df_combined.values*repeatability})

    return df_out


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

    df = df[["mf", "mf_uncertainty"]].sort_index()
    df.index.name = "time"

    return df


def attributes(species, site, network):

    '''
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
        "git_pollution_flag"
        "met_office_baseline_flag"
    '''



def combine_datasets(species, site):

    # Get instructions on how to combine datasets
    with open(paths.root / "data/data_selector.json") as f:
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


