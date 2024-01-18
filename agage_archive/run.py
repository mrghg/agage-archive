import pandas as pd
from shutil import rmtree
from zipfile import ZipFile

from agage_archive import Paths, open_data_file, data_file_list, data_file_path
from agage_archive.data_selection import read_release_schedule, read_data_combination
from agage_archive.io import combine_datasets, combine_baseline, \
    read_nc, read_baseline, read_ale_gage, \
    output_dataset


def run_individual_instrument(network, instrument,
                              verbose = False,
                              baseline = ""):
    """Process individual data files for a given instrument.
    Reads the release schedule for the instrument

    Args:
        instrument (str): Instrument to process. Must match sheet names in release schedule, e.g.:
            "AGAGE", "ALE", "GAGE", "GCMD", ...
        verbose (bool): Print progress to screen
        baseline (str): Baseline flag to use. If empty, don't process baselines
    """

    rs = read_release_schedule(network, instrument)

    if instrument.upper() == "ALE" or instrument.upper() == "GAGE":
        read_function = read_ale_gage
        read_baseline_function = read_baseline
        instrument_out = instrument.upper() + "-GCMD"
    else:
        read_function = read_nc
        read_baseline_function = read_baseline
        instrument_out = instrument.upper()

    # Process for all species and sites
    for species in rs.index:
        for site in rs.columns:
            if rs.loc[species, site].lower() != "x":

                ds = read_function(network, species, site, instrument,
                                verbose=verbose)

                if baseline:
                    ds_baseline = read_baseline_function(network, species, site, instrument,
                                                flag_name = baseline,
                                                verbose = verbose)

                # If multiple instruments, store individual file in subdirectory
                instrument_dates = read_data_combination(network, species, site,
                                                        verbose=False)
                if len(instrument_dates) > 1:
                    output_subpath = f"{species}/individual"
                else:
                    output_subpath = species

                output_dataset(ds, network, instrument=instrument_out,
                               output_subpath=output_subpath,
                               end_date=rs.loc[species, site],
                               verbose=verbose)

                if baseline:
                    if (ds_baseline.time != ds.time).any():
                        raise ValueError(f"Baseline and data files for {species} at {site} have different timestamps")
                    output_dataset(ds_baseline, network, instrument=instrument_out,
                               output_subpath=output_subpath + "/baseline_flags",
                               end_date=rs.loc[species, site],
                               extra="-git-baseline",
                               verbose=verbose)


def run_combined_instruments(network,
                             baseline = False,
                             verbose = False):
    """Process combined data files for a given network.
    Reads the data selection file to determine which sites to process

    Args:
        network (str): Network for output filenames
        baseline (bool): Process baselines. Boolean as only one baseline flag is available (GIT)
        verbose (bool): Print progress to screen
    """

    with open_data_file("data_combination.xlsx", network=network) as data_selection_path:
        sites = pd.ExcelFile(data_selection_path).sheet_names

    for site in sites:

        print(f"Processing files for {site}")

        # Read data selection file
        with open_data_file("data_combination.xlsx", network=network) as f:
            df = pd.read_excel(f,
                            comment="#",
                            sheet_name=site,
                            index_col="Species")

        # Loop through species in index
        for species in df.index:

            # Produce combined dataset
            if verbose:
                print(f"... combining datasets for {species} at {site}")
            ds = combine_datasets(network, species, site, verbose=verbose)

            if baseline:
                if verbose:
                    print(f"... combining baselines for {species} at {site}")
                # Note that GIT baselines is hard-wired here because Met Office not available for ALE/GAGE
                ds_baseline = combine_baseline(network, species, site,
                                            verbose=verbose)

            if verbose:
                print(f"... outputting combined dataset for {species} at {site}")
            output_dataset(ds, network,
                           output_subpath=species,
                           instrument="combined",
                           verbose=verbose)
            
            if baseline:
                if verbose:
                    print(f"... outputting combined baseline for {species} at {site}")
                output_dataset(ds_baseline, network,
                               output_subpath=species + "/baseline_flags",
                               instrument="combined",
                               extra="-git-baseline",
                               verbose=verbose)


def run_all(network,
            delete = True,
            combined = True,
            baseline = True,
            include = [],
            exclude = ["GCPDD"]):
    """Process data files for multiple instruments. Reads the release schedule to determine which
    instruments to process

    Args:
        delete (bool): Delete all files in output directory before running
        combined (bool): Process combined data files
        include (list): List of instruments to process. If empty, process all instruments
        exclude (list): List of instruments to exclude from processing
    """

    if not network:
        raise ValueError("Must specify network")

    path = Paths(network, errors="ignore")

    out_pth = data_file_path("", network=network, sub_path=path.output_path, errors="ignore")

    if delete:
        # Clear output directory, removing all files and subdirectories
        network, sub_path, files = data_file_list(network=network,
                                                  sub_path=path.output_path,
                                                  errors="ignore")
        
        print(f'Deleting all files in {out_pth}')

        if out_pth.suffix == ".zip" and out_pth.exists():

            out_pth.unlink()

        else:

            for f in files:

                pth = out_pth / f

                # Make sure pth is in data/network directory (for safety)
                if pth.parents[2] == path.data and pth.parents[1].name == network:
                    if pth.is_file():
                        pth.unlink()
                    elif pth.is_dir():
                        rmtree(pth)
                else:
                    print(f"Warning: {pth} must be in a data/network directory")

    # If out_pth is a zip file that doesn't exist, create it
    if out_pth.suffix == ".zip" and not out_pth.exists():
        with ZipFile(out_pth, "w") as f:
            pass

    # Must run combined instruments first
    if combined:
        run_combined_instruments(network, baseline=baseline, verbose=True)

    # If include is empty, process all instruments in release schedule
    if len(include) == 0:
        with open_data_file("data_release_schedule.xlsx", network=network) as frs:
            instruments = pd.ExcelFile(frs).sheet_names
    else:
        instruments = include

    # Processing
    for instrument in instruments:
        if instrument not in exclude:
            baseline_flag = {True: "git_pollution_flag", False: ""}[baseline]
            run_individual_instrument(network, instrument, 
                                    baseline=baseline_flag, verbose=True)


if __name__ == "__main__":

    run_all("agage")