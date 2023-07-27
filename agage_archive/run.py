import pandas as pd

from agage_archive import Paths
from agage_archive.data_selection import read_release_schedule
from agage_archive.io import combine_datasets, read_agage, read_ale_gage, output_dataset

path = Paths()


def run_individual_instrument(instrument,
                              verbose=False):
    """Process individual data files for a given instrument.
    Reads the release schedule for the instrument

    Args:
        instrument (str): Instrument to process. Must match sheet names in release schedule, e.g.:
            "AGAGE", "ALE", "GAGE", "GCMD", ...
        verbose (bool): Print progress to screen
    """

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
                

def run_combined_instruments(network = "AGAGE",
                             verbose = False):
    """Process combined data files for a given network.
    Reads the data selection file to determine which sites to process

    Args:
        network (str): Network for output filenames
        verbose (bool): Print progress to screen
    """

    file_path = path.root / "data" / "data_selection" / "data_selection.xlsx"

    # Read sheet names in file_path to determine which sites to process
    sites = pd.ExcelFile(file_path).sheet_names

    for site in sites:

        print(f"Processing files for {site}")

        df = pd.read_excel(file_path,
                        comment="#",
                        sheet_name=site,
                        index_col="Species")

        # Loop through species in index
        for species in df.index:

            # Produce combined dataset
            if verbose:
                print(f"... combining datasets for {species} at {site}")
            ds = combine_datasets(species, site, verbose=verbose)

            if verbose:
                print(f"... outputting combined dataset for {species} at {site}")
            output_dataset(ds, network, instrument="combined", verbose=verbose)