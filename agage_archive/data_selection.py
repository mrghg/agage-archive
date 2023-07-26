import pandas as pd
import numpy as np

from agage_archive import Paths
from agage_archive.formatting import format_species


paths = Paths()


def read_instrument_dates_xlsx(species, site):
    '''Read instrument dates from Excel file

    Args:
        species (str): Species
        site (str): Site code

    Returns:
        dict: Dictionary of instrument dates
    '''

    path = paths.root / "data" / "data_selection" / "data_selection.xlsx"

    warning_message = f"WARNING: No instrument dates found for {species} at {site.upper()}. Assuming GCMS-Medusa"

    default_output = {"GCMS-Medusa": [None, None]}

    # Determine if data_exclude.xlsx contains sheet name called site.upper()
    if site.upper() not in pd.ExcelFile(path).sheet_names:
        print(warning_message)
        return default_output

    # Read data_selection
    df = pd.read_excel(path,
                    comment="#",
                    sheet_name=site.upper(),
                    index_col="Species")
    
    # Look for species name in table, return Medusa if not there
    df = df[df.index == format_species(species)]

    if len(df) == 0:
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


def read_release_schedule(instrument,
                          species = None,
                          site = None):

    with open(paths.root / "data/data_selection/data_release_schedule.xlsx", "+rb") as f:
        df_all = pd.read_excel(f, sheet_name=instrument)

        # Get index of row that contains "General release date"
        idx = df_all[df_all.iloc[:, 0] == "General release date"].index[0]
        general_end_date = df_all.iloc[idx, 1]

        df = pd.read_excel(f, sheet_name=instrument, skiprows=idx+2)

    # Remove whitespace
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # If NaN in df, replace with general end date
    df = df.fillna(general_end_date)

    df.set_index("Species", inplace=True)

    if species is not None:
        return df.loc[species, site]
    else:
        return df


def calibration_scale_default(species):
    '''Get default calibration scale

    Args:
        species (str): Species

    Returns:
        str: Calibration scale
    '''

    # Read scale_defaults csv file
    scale_defaults = pd.read_csv(paths.root / "data/scale_defaults.csv",
                                index_col="Species")
    
    scale_defaults = scale_defaults.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    if format_species(species) not in scale_defaults.index.str.lower():
        return scale_defaults.loc["all", "calibration_scale"]
    else:
        return scale_defaults.loc[format_species(species), "calibration_scale"]


def data_exclude(ds, species, site, instrument):
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
    if site.upper() not in pd.ExcelFile(paths.root / "data/data_selection/data_exclude.xlsx").sheet_names:
        return ds

    # Read data_exclude
    data_exclude = pd.read_excel(paths.root / "data/data_selection/data_exclude.xlsx",
                                comment="#",
                                sheet_name=site.upper())
    
    # Remove whitespace from strings
    data_exclude = data_exclude.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # columns are Species, Instrument, Start, End
    # Find rows that match species and instrument and then extract start and end dates
    data_exclude = data_exclude[(data_exclude["Species"] == format_species(species)) &
                                (data_exclude["Instrument"] == instrument)][["Start", "End"]]

    if len(data_exclude) == 0:
        return ds
    else:
        # Replace values in xarray dataset with NaN between start and end dates
        for start, end in data_exclude.values:
            ds.loc[dict(time=slice(start, end))] = np.nan
        return ds