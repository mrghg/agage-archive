import pandas as pd
from shutil import rmtree

from agage_archive import Paths
from agage_archive.data_selection import read_release_schedule, read_data_combination
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

                # Determine which read function to use
                if instrument == "ALE" or instrument == "GAGE":
                    ds = read_function(species, site, instrument,
                                       verbose=verbose)
                else:
                    ds = read_function(species, site, instrument,
                                       verbose=verbose)

                # If multiple instruments, store individual file in subdirectory
                instrument_dates = read_data_combination(species, site, verbose=False)
                if len(instrument_dates) > 1:
                    output_subpath = f"{species}/individual"
                else:
                    output_subpath = species

                output_dataset(ds, network, instrument=instrument_out,
                               output_subpath=output_subpath,
                               end_date=rs.loc[species, site],
                               verbose=verbose)
                

def run_combined_instruments(network = "AGAGE",
                             verbose = False):
    """Process combined data files for a given network.
    Reads the data selection file to determine which sites to process

    Args:
        network (str): Network for output filenames
        verbose (bool): Print progress to screen
    """

    data_selection_path = path.root / "data" / "data_selection" / "data_combination.xlsx"

    # Read sheet names in file_path to determine which sites to process
    sites = pd.ExcelFile(data_selection_path).sheet_names

    for site in sites:

        print(f"Processing files for {site}")

        df = pd.read_excel(data_selection_path,
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
            output_dataset(ds, network,
                           output_subpath=species,
                           instrument="combined",
                           verbose=verbose)


def run_all(delete = True):

    #TODO: SORT THIS OUT!!!!!!
    exclude = ["GCPDD"]

    rs_path = path.root / "data" / "data_selection" / "data_release_schedule.xlsx"

    if delete:
        print("Deleting all files in output directory")
        # Clear output directory, removing all files and subdirectories
        pths = path.output.glob("*")

        for pth in pths:
            if pth.is_file():
                pth.unlink()
            elif pth.is_dir():
                rmtree(pth)

    # Must run combined instruments first
    run_combined_instruments(verbose=True)

    instruments = pd.ExcelFile(rs_path).sheet_names

    for instrument in instruments:

        if instrument not in exclude:
            run_individual_instrument(instrument, verbose=True)


if __name__ == "__main__":

    run_all()