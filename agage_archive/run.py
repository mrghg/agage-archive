import pandas as pd

from agage_archive import Paths
from agage_archive.io import combine_datasets, read_agage, read_ale_gage, output_dataset

path = Paths()


def read_release_schedule(instrument):

    with open(path.root / "data/data_selection/data_release_schedule.xlsx", "+rb") as f:
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

    return df


def run_individual_instrument(instrument):

    rs = read_release_schedule(instrument)

    if instrument == "ALE" or instrument == "GAGE":
        network = instrument
        instrument_out = "GCMD"
        read_function = read_ale_gage
    else:
        network = "AGAGE"
        instrument_out = instrument
        read_function = read_agage

    # Process for all species and sites
    for species in rs.index:
        for site in rs.columns:
            if rs.loc[species, site].lower() != "x":
                if instrument == "ALE" or instrument == "GAGE":
                    ds = read_function(species, site, instrument)
                else:
                    ds = read_function(species, site, instrument)
                output_dataset(ds, network, instrument=instrument_out,
                               end_date=rs.loc[species, site])