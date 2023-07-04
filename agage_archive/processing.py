import json
import pandas as pd
import xarray as xr
import numpy as np
from datetime import datetime

from agage_archive.io import Paths, read_ale_gage, read_agage
from agage_archive.util import is_number


instrument_number = {"ALE": 0,
                    "GAGE": 1,
                    "GCMD": 2,
                    "GCMS-ADS": 3,
                    "GCMS-Medusa": 4}

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

        instrument_rec.append(instrument)

        if instrument in ["ALE", "GAGE"]:
            ds = read_ale_gage(species, site, instrument)
        else:
            ds = read_agage(species, site, instrument)

        attrs.append(ds.attrs)

        #comments.append(ds.attrs["comment"])

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

    # # Add details on instruments to global attributes
    ds_combined = attributes_instruments(ds_combined,
                                         instrument_rec,
                                         attrs,
                                         dates_rec)

    # # Extend comment attribute describing all datasets
    # ds_combined.attrs["comment"] = "Combined AGAGE/GAGE/ALE dataset combined from the following individual sources: ---- " + \
    #     "|---| ".join(comments)
    
    # Add site code
    ds_combined.attrs["site_code"] = site.upper()

    return ds_combined


def attributes_instruments(ds,
                           instrument_list,
                           attribute_list,
                           dates_list):

    ds.attrs = {}

    # Loop through instruments in reverse
    for instrument in list(instrument_number.keys())[::-1]:
        if len(ds.attrs) == 0:
            if instrument in instrument_list:
                ds.attrs = attribute_list[instrument_list.index(instrument)]
                if not "instrument" in ds.attrs:
                    ds.attrs["instrument"] = instrument
                    ds.attrs["instrument_date"] = dates_list[instrument_list.index(instrument)]
                    ds.attrs["instrument_comment"] = ""
        else:
            if instrument in instrument_list:
                # Find maximum existing instrument number
                instrument_max = 0
                for attr in ds.attrs:
                    if "instrument" in attr:
                        if is_number(attr.split("_")[-1]):
                            if int(attr.split("_")[-1]) > instrument_max:
                                instrument_max = int(attr.split("_")[-1])

                # Add new instrument
                ds.attrs[f"instrument_{instrument_max + 1}"] = instrument
                ds.attrs[f"instrument_{instrument_max + 1}_date"] = dates_list[instrument_list.index(instrument)]
                ds.attrs[f"instrument_{instrument_max + 1}_comment"] = ""

    return ds
