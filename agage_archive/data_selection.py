import pandas as pd
import numpy as np
import json

from agage_archive.config import open_data_file
from agage_archive.formatting import format_species


def read_data_combination(network, species, site,
                        verbose=True):
    '''Read instrument dates from Excel file

    Args:
        species (str): Species
        site (str): Site code
        verbose (bool): Print warning message if no instrument dates found

    Returns:
        dict: Dictionary of instrument dates
    '''

    warning_message = f"WARNING: No instrument dates found for {species} at {site.upper()}. Assuming GCMS-Medusa"

    default_output = {"GCMS-Medusa": [None, None]}

    # Determine if data_exclude.xlsx contains sheet name called site.upper()
    with open_data_file("data_combination.xlsx", network=network) as f:
        if site.upper() not in pd.ExcelFile(f).sheet_names:
            if verbose:
                print(warning_message)
            return default_output

    # Read data_selection
    with open_data_file("data_combination.xlsx", network=network) as f:
        df = pd.read_excel(f,
                        comment="#",
                        sheet_name=site.upper(),
                        index_col="Species")
    
    # Look for species name in table, return Medusa if not there
    df = df[df.index == format_species(species)]

    if len(df) == 0:
        if verbose:
            print(warning_message)
        return default_output
    
    if len(df) > 1: 
        raise ValueError("Can only have each species appear at most once in data_selection table")

    # Extract instrument names from column names
    instruments = [col.split(" ")[0] for col in df.columns]

    # For each instrument, find start and end dates
    instrument_dates = {}

    for instrument in set(instruments):

        dates = df.loc[:, [f"{instrument} start", f"{instrument} end"]].values[0].astype(str)

        if "x" not in dates:
            instrument_dates[instrument] = [date if date != "nan" else None for date in dates]

    return instrument_dates


def read_release_schedule(network, instrument,
                          species = None,
                          site = None,
                          public=True):
    '''Read release schedule from Excel file

    Args:
        network (str): Network
        instrument (str): Instrument
        species (str): Species
        site (str): Site code
        public (bool, optional): Whether the dataset is for public release. Default to True

    Returns:
        pd.DataFrame: Release schedule
    '''

    with open_data_file("data_release_schedule.xlsx", network=network) as f:
        df_all = pd.read_excel(f, sheet_name=instrument)

        # Get index of row that contains "General release date"
        idx = df_all[df_all.iloc[:, 0] == "General release date"].index[0]
        general_end_date = df_all.iloc[idx, 1]

        # Determine if general_end_date is a string
        if not isinstance(general_end_date, str):
            print(f"WARNING: No general end date found for {instrument}. Assuming no limit.")

        df = pd.read_excel(f, sheet_name=instrument, skiprows=idx+2)

    # Remove whitespace
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)

    # If NaN in df, replace with general end date
    df = df.fillna(general_end_date)
    
    # For private release data set general end date far into future
    if not public:
        df = df.where(df.values != general_end_date, "2100-01-01 00:00")

    df.set_index("Species", inplace=True)

    if species is not None:

        if format_species(species) not in df.index.str.lower():
            raise ValueError(f"No release schedule found for {species} at {site}")

        if df.loc[species, site] == "x":
            # Return a value before any data was collected, to remove everything
            return "1970-01-01"
        elif not isinstance(df.loc[species, site], str):
            return None
        else:
            return df.loc[species, site]
    else:
        return df


def calibration_scale_default(network, species):
    '''Get default calibration scale

    Args:
        species (str): Species

    Returns:
        str: Calibration scale
    '''

    # Read scale_defaults csv file
    with open_data_file("scale_defaults.csv", network=network) as f:
        scale_defaults = pd.read_csv(f, index_col="Species")
    
    scale_defaults = scale_defaults.map(lambda x: x.strip() if isinstance(x, str) else x)

    if format_species(species) not in scale_defaults.index.str.lower():
        raise ValueError(f"No default calibration scale found for {species}")
    else:
        return scale_defaults.loc[format_species(species), "calibration_scale"]


def read_data_exclude(ds, species, site, instrument):
    '''Read data_exclude file and return start and end date for exclusion

    Args:
        ds (xr.Dataset): Dataset
        species (str): Species
        site (str): Site code
        instrument (str): Instrument

    Returns:
        xr.Dataset: Dataset with NaNs between start and end dates
    '''

    # Determine if data_exclude.xlsx contains sheet name called site.upper()
    with open_data_file("data_exclude.xlsx", network = ds.attrs["network"]) as f:
        if site.upper() not in pd.ExcelFile(f).sheet_names:
            return ds

    # Read data_exclude
    with open_data_file("data_exclude.xlsx", network = ds.attrs["network"]) as f:
        data_exclude = pd.read_excel(f,
                                    comment="#",
                                    sheet_name=site.upper())
    
    # Read variable defaults to find what to do with missing data
    with open_data_file("variables.json", this_repo=True) as f:
        variable_defaults = json.load(f)

    # Read variable defaults for non-public data to find what to do with missing data
    with open_data_file("variables_not_public.json", this_repo=True) as f:
        variable_np_defaults = json.load(f)
    variable_defaults.update(variable_np_defaults)

    # Remove whitespace from strings
    data_exclude = data_exclude.map(lambda x: x.strip() if isinstance(x, str) else x)

    # columns are Species, Instrument, Start, End
    # Find rows that match species and instrument and then extract start and end dates
    data_exclude = data_exclude[(data_exclude["Species"] == format_species(species)) &
                                (data_exclude["Instrument"] == instrument)][["Start", "End"]]

    if len(data_exclude) == 0:
        return ds
    else:
        # Replace values in xarray dataset with NaN between start and end dates
        # The remove_flagged value in variables.json determines what to do with the missing data
        for start, end in data_exclude.values:
            exclude_dict = dict(time = slice(start, end))
            for var in ds.variables:
                if "time" in var:
                    continue
                elif variable_defaults[var]["remove_flagged"] == "True":
                    ds[var].loc[exclude_dict] = np.nan
                elif variable_defaults[var]["remove_flagged"] == "False":
                    continue
                elif variable_defaults[var]["remove_flagged"] == "Zero":
                    ds[var].loc[exclude_dict] = 0
                else:
                    raise ValueError(f"Unknown remove_flagged value for {var}. Check in variables.json or variables_not_public.json")

        return ds