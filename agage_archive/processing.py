import json
import pandas as pd
import xarray as xr
import numpy as np
from datetime import datetime
import os
import configparser

from agage_archive import Paths, get_path
from agage_archive.util import is_number


nc4_types = {"f4": "float32",
            "f8": "float64",
            "i4": "int32",
            "i2": "int16",
            "i8": "int64",
            "i1": "int8",
            "u1": "uint8",
            "u2": "uint16",
            "u4": "uint32",
            "u8": "uint64",
            "S1": "S1"}


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


def format_attributes_global_instrument(ds,
                                instrument,
                                return_attributes=False):
    '''Add global attributes for instrument

    Args:
        ds (xr.Dataset): Dataset
        instrument (str): Instrument
    '''

    attrs = ds.attrs.copy()

    if "instrument" not in ds.attrs.keys():
        attrs["instrument"] = instrument
    if "instrument_date" not in ds.attrs.keys():
        attrs["instrument_date"] = f"{ds.time[0].dt.strftime('%Y-%m-%d').values}"
    if "instrument_comment" not in ds.attrs.keys():
        attrs["instrument_comment"] = ""

    if return_attributes:
        return {key: attrs[key] for key in attrs if "instrument" in key}
    else:
        ds_out = ds.copy(deep=True)
        ds_out.attrs = attrs.copy()
        return ds_out


def format_attributes_global_instruments(ds,
                                        instruments = [],
                                        return_attributes=False):
    '''Combine instrument details in global attributes

    Args:
        ds (xr.Dataset): Dataset
        instruments (list): List of instrument dictionaries containing the keys "instrument", "instrument_date" and "instrument_comment"
        
    Returns:
        xr.Dataset: Dataset with updated global attributes
    '''

    # If no instruments specified, check if they are in the dataset attributes
    if len(instruments) == 0:
        if "instrument" in ds.attrs:
            instruments = [{key: ds.attrs[key] for key in ds.attrs if "instrument" in key}]
        else:
            raise ValueError("No instrument* attributes specified and none found in dataset attributes. " + \
                             "As a minimum, set keyword instrument = [{'instrument': 'INSTRUMENT_NAME'}]")

    # If only one instrument, make sure formatted properly and return
    elif len(instruments) == 1:
        return format_attributes_global_instrument(ds, instruments[0]["instrument"],
                                                return_attributes=return_attributes)

    attrs = {}

    # Remove instrument attributes so that we can repopulate them
    for attr in ds.attrs:
        if "instrument" not in attr:
            attrs[attr] = ds.attrs[attr]

    dates = []
    suffixes = []

    for instrument in instruments:
        suffix = ["",]

        has_instrument = False
        has_instrument_date = False

        for key, value in instrument.items():
            # For now, just get the date of the first instrument, as it makes sorting easier
            # There can be multiple, but it's probably not a big deal if they aren't in exactly the right order
            if key == "instrument_date":
                dates.append(value)
                has_instrument_date = True
            
            if is_number(key.split("instrument_")[-1]):
                # Store numerical suffix, prepending with underscore
                suffix.append("_" + key.split("instrument_")[-1])
                has_instrument = True

            if key == "instrument":
                has_instrument = True
                
        suffixes.append(suffix)

        if not has_instrument:
            raise ValueError("No instrument attribute found")
        if not has_instrument_date:
            raise ValueError("No instrument_date attribute found")
        
    #Sort index by date
    idx = np.argsort(dates)[::-1]
    # Apply this sort order to instruments dictionary
    instruments_sorted = [instruments[i] for i in idx]
    suffixes_sorted = [suffixes[i] for i in idx]

    instrument_count = 0

    # Loop through instruments in reverse
    for instrument, suffix in zip(instruments_sorted, suffixes_sorted):

        for n in suffix:

            # Relabel instrument number
            if instrument_count == 0:
                suffix_new = ""
            else:
                suffix_new = "_" + str(instrument_count)

            for attr in ["instrument", "instrument_date", "instrument_comment"]:
                if attr + n in instrument.keys():
                    attrs[attr + suffix_new] = instrument[attr + n]
                else:
                    print("WARNING: No " + attr + " found for instrument " + n + ". Setting to empty string")
                    attrs[attr + suffix_new] = ""
        
            instrument_count += 1

    if return_attributes:
        return {key: attrs[key] for key in attrs if "instrument" in key}
    else:
        ds_out = ds.copy(deep=True)
        ds_out.attrs = attrs.copy()
        return ds_out


def format_dataset(ds,
                variable_translate={},
                instruments = [],
                species = None,
                units = None,
                calibration_scale = None):
    '''Format attributes, variables and encoding

    Args:
        ds (xr.Dataset): Dataset
    
    Returns:
        xr.Dataset: Dataset with updated global attributes
    '''

    ds = format_attributes(ds,
                        instruments = instruments,
                        species = species,
                        units = units,
                        calibration_scale = calibration_scale)
    
    ds = format_variables(ds,
                        variable_translate = variable_translate,
                        species = species,
                        units = units,
                        calibration_scale = calibration_scale)

    #TODO: Format comment string?

    return ds


def format_variables(ds,
                    variable_translate={},
                    species = None,
                    units = None,
                    calibration_scale = None):
    '''Standardise variable names and units

    Args:
        ds (xr.Dataset): Dataset
        variable_translate (dict, optional): Dictionary of variable translations. Defaults to {}.
        species (str, optional): Species name. Defaults to None, in which case it is looked up in the dataset attributes.
        units (str, optional): Units. Defaults to None, in which case it is looked up in the dataset attributes.
        calibration_scale (str, optional): Calibration scale. Defaults to None, in which case it is looked up in the dataset attributes.

    Returns:
        xr.Dataset: Dataset with standardised variable names and units

    Raises:
        ValueError: If variable not found in dataset

    '''

    with open(paths.root / "data/variables.json") as f:
        variables = json.load(f)

    attrs = ds.attrs.copy()
    
    vars_out = {}

    # Loop through standard variable names
    for var in variables:

        # Do we need to translate any variable names?
        if var in variable_translate:
            var_ds = variable_translate[var]
        else:
            var_ds = var
        
        # Copy and convert variable from dataset
        if var_ds in ds.variables:
            # Convert type to that contained in variables.json
            typ = np.__getattribute__(nc4_types[variables[var]["encoding"]["dtype"]])
            
            vars_out[var] = ("time", ds[var_ds].values.copy().astype(typ))

        else:
            if variables[var]["optional"] == "False":
                raise ValueError(f"Variable {var_ds} not found in dataset. " + \
                                "Use variable_translate to map to a different variable name.")
    
    # Time can't be in variable list
    del vars_out["time"]

    # Create new dataset
    ds = xr.Dataset(vars_out,
                    coords = {"time": ds.time.copy()}, 
                    attrs=attrs)

    # Add variable attributes and encoding
    for var in ds.variables:
        ds[var].attrs = variables[var]["attrs"]
        ds[var].encoding = variables[var]["encoding"]

        # for mole fractions, replace % with species name, etc.
        if "mf" in var:
            ds[var].attrs["long_name"] = ds[var].attrs["long_name"].replace("%",
                                                        lookup_locals_and_attrs("species", locals(), attrs))
            ds[var].attrs["standard_name"] = ds[var].attrs["standard_name"].replace("%",
                                                        lookup_locals_and_attrs("species", locals(), attrs))
            ds[var].attrs["units"] = lookup_locals_and_attrs("units", locals(), attrs)
            ds[var].attrs["calibration_scale"] = lookup_locals_and_attrs("calibration_scale", locals(), attrs)

    return ds


def lookup_locals_and_attrs(v, local, attrs):
    '''Look up variable in locals and attrs, and format the output

    Args:
        v (str): Variable name
        local (dict): Dictionary of local variables
        attrs (dict): Dictionary of attributes

    Returns:
        str: Variable name
    '''

    #TODO: this assumes that the variable is in locals
    if local[v] is None:
        # If not set as keyword, check if it is in dataset attributes
        if v not in attrs:
            raise ValueError(f"{v} not set and not found in dataset attributes")
        else:
            # If in dataset attributes, format it
            return eval(f"format_{v}('{attrs[v]}')")
    else:
        # If set, format it
        return eval(f"format_{v}('{local[v]}')")


def format_attributes(ds, instruments = [],
                      species = None,
                      units = None,
                      calibration_scale = None):
    '''Format attributes

    Args:
        ds (xr.Dataset): Dataset
        instruments (list, optional): List of instrument dictionaries containing the keys "instrument", "instrument_date" and "instrument_comment"

    Returns:
        xr.Dataset: Dataset with formatted attributes
    '''

    with open(paths.root / "data/attributes.json") as f:
        attributes_default = json.load(f)

    attrs = {}

    for attr in attributes_default:

        if "instrument" in attr:
            # Update instrument attributes
            attrs.update(format_attributes_global_instruments(ds,
                                                        instruments,
                                                        return_attributes=True))

        elif attr == "file_created_by":
            # Get username
            attrs[attr] = lookup_username()

        elif attr == "file_created":
            # Get current time
            attrs[attr] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        else:
            # Set default
            attrs[attr] = attributes_default[attr]

            # If attribute is in dataset, and not empty, overwrite
            if attr in ds.attrs:
                if ds.attrs[attr] != "":
                    attrs[attr] = ds.attrs[attr]

    # Format certain key attributes, and determine if they have been set as keywords
    for v in ["species", "units", "calibration_scale"]:
        attrs[v] = lookup_locals_and_attrs(v, locals(), ds.attrs.copy())

    ds_out = ds.copy(deep=True)
    ds_out.attrs = attrs.copy()

    return ds_out


def format_species(species):
    '''Format species name

    Args:
        species (str): Species name

    Returns:
        str: Formatted species name
    '''

    return species.lower() #.replace("-", "")


def format_units(units):
    '''Format units

    Args:
        units (str): Units

    Returns:
        str: Formatted units
    '''

    unit_translator = {"ppm": "1e-6",
                       "ppb": "1e-9",
                       "ppt": "1e-12",
                       "ppq": "1e-15",
                       "nmol/mol": "1e-9",
                       "nmol mol-1": "1e-9",
                       "pmol/mol": "1e-12",
                       "pmol mol-1": "1e-12",
                       }

    if units in unit_translator:
        return unit_translator[units]
    else:
        return units


def format_calibration_scale(scale):
    '''Format scale

    Args:
        scale (str): Scale

    Returns:
        str: Formatted scale
    '''

    return scale


def lookup_username():
    '''Look up username

    Returns:
        str: Username
    '''

    # Check if config file exists
    config_file = get_path("config.ini")
    if not config_file.exists():
        raise FileNotFoundError(
            "Config file not found. Try running util.setup first")
    
    # Take username from config file if it exists, otherwise try to get it from the system
    config = configparser.ConfigParser()
    config.read(config_file)

    if "User" in config:
        return config["User"]["name"]
    else:
        try:
            return os.environ["USER"]
        except:
            try:
                return os.environ["USERNAME"]
            except:
                try:
                    return os.environ["LOGNAME"]
                except:
                    return "unknown user"