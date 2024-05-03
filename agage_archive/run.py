import pandas as pd
from shutil import rmtree
from zipfile import ZipFile

from agage_archive.config import Paths, open_data_file, data_file_list, data_file_path, copy_to_archive
from agage_archive.data_selection import read_release_schedule, read_data_combination
from agage_archive.io import combine_datasets, combine_baseline, \
    read_nc, read_baseline, read_ale_gage, read_gcwerks_flask, \
    output_dataset
from agage_archive.convert import monthly_baseline


def run_timestamp_checks(ds,
                        ds_baseline=None,
                        species="",
                        site=""):

    # Check for duplicate time stamps
    if ds["time"].to_series().duplicated().any():
        raise ValueError(f"Duplicate timestamps in {species} at {site}")
    if ds_baseline:
        if ds_baseline["time"].to_series().duplicated().any():
            raise ValueError(f"Duplicate timestamps in baseline for {species} at {site}")
    
    # check that the time stamps are the same in the data and baseline files
    if ds_baseline:
        if (ds_baseline.time != ds.time).any():
            raise ValueError(f"Data and baseline files for {species} at {site} have different timestamps")
        

def run_individual_instrument(network, instrument,
                              verbose = False,
                              baseline = "",
                              monthly = False,
                              species = [],
                              public=True):
    """Process individual data files for a given instrument.
    Reads the release schedule for the instrument

    Args:
        instrument (str): Instrument to process. Must match sheet names in release schedule, e.g.:
            "AGAGE", "ALE", "GAGE", "GCMD", ...
        verbose (bool): Print progress to screen
        baseline (str): Baseline flag to use. If empty, don't process baselines
        monthly (bool): Produce monthly baseline files
        species (list): List of species to process. If empty, process all species
        public (bool, optional): Whether the dataset is for public release. Default to True.
    """
    
    rs = read_release_schedule(network, instrument, public=public)

    if instrument.upper() == "ALE" or instrument.upper() == "GAGE":
        read_function = read_ale_gage
        read_baseline_function = read_baseline
        instrument_out = instrument.upper() + "-GCMD"
    elif instrument.upper() == "GCMS-MEDUSA-FLASK":
        read_function = read_gcwerks_flask
        instrument_out = "GCMS-MEDUSA-FLASK"
    else:
        read_function = read_nc
        read_baseline_function = read_baseline
        instrument_out = instrument.upper()

    if species:
        # Process only those species that are in the release schedule
        species_to_process = [sp for sp in species if sp in rs.index.values]
        if not species_to_process:
            print(f"No species to process for {instrument}, skipping...")
            return
    else:
        # Process all species in the release schedule
        species_to_process = rs.index.values

    error_log = []

    # Process for all species and sites
    for sp in species_to_process:
        for site in rs.columns:

            try:

                if rs.loc[sp, site].lower() != "x":

                    ds = read_function(network, sp, site, instrument,
                                    public = public, verbose=verbose)

                    if baseline:
                        ds_baseline = read_baseline_function(network, sp, site, instrument,
                                                    flag_name = baseline,
                                                    verbose = verbose)
                    else:
                        ds_baseline = None
                        
                    run_timestamp_checks(ds, ds_baseline, sp, site)

                    # If multiple instruments, store individual file in subdirectory
                    instrument_dates = read_data_combination(network, sp, site,
                                                            verbose=False)
                    if len(instrument_dates) > 1:
                        output_subpath = f"event/{sp}/individual"
                    else:
                        output_subpath = f"event/{sp}"

                    output_dataset(ds, network, instrument=instrument_out,
                                output_subpath=output_subpath,
                                end_date=rs.loc[sp, site],
                                public=public,
                                verbose=verbose)

                    if baseline:
                        if (ds_baseline.time != ds.time).any():
                            raise ValueError(f"Baseline and data files for {sp} at {site} have different timestamps")
                        output_dataset(ds_baseline, network, instrument=instrument_out,
                                output_subpath=output_subpath + "/baseline_flags",
                                end_date=rs.loc[sp, site],
                                extra="-git-baseline",
                                public=public,
                                verbose=verbose)
                        
                        if monthly:
                            ds_baseline_monthly = monthly_baseline(ds, ds_baseline)
                            output_dataset(ds_baseline_monthly, network, instrument=instrument_out,
                                output_subpath=output_subpath.replace("event", "monthly"),
                                end_date=rs.loc[sp, site],
                                extra="-monthly",
                                public=public,
                                verbose=verbose)

                    else:
                        if monthly:
                            raise NotImplementedError("Monthly baseline files can only be produced if baseline flag is specified")

            except Exception as e:
                error_log.append((site, sp, e))

    if error_log:
        # save errors to file
        with open(data_file_path("error_log_individual.txt", network=network, errors="ignore"), "w") as f:
            # write the date and time of the error
            f.write("Processing attempted on " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
            for error in error_log:
                f.write(f"{error[0]} {error[1]}: {error[2]}\n")
        
        print("!!! Errors occurred during processing. See error_log_individual.txt for details")


def run_combined_instruments(network,
                             baseline = False,
                             monthly = False,
                             verbose = False,
                             species = [],
                             public = True):
    """Process combined data files for a given network.
    Reads the data selection file to determine which sites to process

    Args:
        network (str): Network for output filenames
        baseline (bool): Process baselines. Boolean as only one baseline flag is available (GIT)
        monthly (bool): Produce monthly baseline files
        verbose (bool): Print progress to screen
        species (list): List of species to process. If empty, process all species
        public (bool, optional): Whether the dataset is for public release. Default to True.
    """

    with open_data_file("data_combination.xlsx", network=network) as data_selection_path:
        sites = pd.ExcelFile(data_selection_path).sheet_names

    # Create error log
    error_log = []

    for site in sites:

        print(f"Processing files for {site}")

        # Read data selection file
        with open_data_file("data_combination.xlsx", network=network) as f:
            df = pd.read_excel(f,
                            comment="#",
                            sheet_name=site,
                            index_col="Species")

        # Determine species to process
        if species:
            # Process only those species that are in the data selection file
            species_to_process = [sp for sp in species if sp in df.index.values]
            if not species_to_process:
                print(f"No species to process for {site}, skipping...")
                continue
        else:
            # Process all species in the data selection file
            species_to_process = df.index.values

        # Loop through species in index
        for species in species_to_process:

            try:

                # Produce combined dataset
                if verbose:
                    print(f"... combining datasets for {species} at {site}")
                ds = combine_datasets(network, species, site,
                                    verbose=verbose, public=public)

                if baseline:
                    if verbose:
                        print(f"... combining baselines for {species} at {site}")
                    # Note that GIT baselines is hard-wired here because Met Office not available for ALE/GAGE
                    ds_baseline = combine_baseline(network, species, site,
                                                verbose=verbose, public=public)

                else:
                    ds_baseline = None

                # Check for duplicate time stamps
                run_timestamp_checks(ds, ds_baseline, species, site)

                output_subpath = f"event/{species}"

                if verbose:
                    print(f"... outputting combined dataset for {species} at {site}")
                output_dataset(ds, network,
                            output_subpath=output_subpath,
                            instrument="combined",
                            public=public,
                            verbose=verbose)
                
                if baseline:
                    if verbose:
                        print(f"... outputting combined baseline for {species} at {site}")
                    output_dataset(ds_baseline, network,
                                output_subpath=output_subpath + "/baseline_flags",
                                instrument="combined",
                                extra="-git-baseline",
                                public=public,
                                verbose=verbose)

                    if monthly:
                        ds_baseline_monthly = monthly_baseline(ds, ds_baseline)
                        output_dataset(ds_baseline_monthly, network,
                                output_subpath=output_subpath.replace("event", "monthly"),
                                instrument="combined",
                                extra="-monthly",
                                public=public,
                                verbose=verbose)

                else:
                    if monthly:
                        raise NotImplementedError("Monthly baseline files can only be produced if baseline flag is specified")

            except Exception as e:
                error_log.append((site, species, e))
    
    if error_log:
        # save errors to file
        with open(data_file_path("error_log_combined.txt", network=network, errors="ignore"), "w") as f:
            # write the date and time of the error
            f.write("Processing attempted on " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
            for error in error_log:
                f.write(f"{error[0]} {error[1]}: {error[2]}\n")
        
        print("!!! Errors occurred during processing. See error_log_combined.txt for details")


def run_all(network,
            delete = True,
            combined = True,
            baseline = True,
            monthly = True,
            instrument_include = [],
            instrument_exclude = ["GCPDD"],
            species = [],
            public = True):
    """Process data files for multiple instruments. Reads the release schedule to determine which
    instruments to process

    Args:
        delete (bool): Delete all files in output directory before running
        combined (bool): Process combined data files
        include (list): List of instruments to process. If empty, process all instruments
        exclude (list): List of instruments to exclude from processing
        baseline (bool): Process baselines. Boolean as only one baseline flag is available (GIT)
        monthly (bool): Produce monthly baseline files
        verbose (bool): Print progress to screen
        species (list): List of species to process. If empty, process all species
        public (bool, optional): Whether the dataset is for public release. Default to True.
        
    """

    if not network:
        raise ValueError("Must specify network")

    if not isinstance(network, str):
        raise TypeError("network must be a string")
    
    if not isinstance(delete, bool):
        raise TypeError("delete must be a boolean")
    
    if not isinstance(combined, bool):
        raise TypeError("combined must be a boolean")
    
    if not isinstance(baseline, bool):
        raise TypeError("baseline must be a boolean")
    
    if not isinstance(monthly, bool):
        raise TypeError("monthly must be a boolean")
    
    if not isinstance(instrument_include, list):
        raise TypeError("instrument_include must be a list")
    
    if not isinstance(instrument_exclude, list):
        raise TypeError("instrument_exclude must be a list")
    
    if not isinstance(species, list):
        raise TypeError("species must be a list")

    path = Paths(network, errors="ignore")

    # Delete log files, if they exist
    for log_file in ["error_log_combined.txt", "error_log_individual.txt"]:
        try:
            data_file_path(log_file, network=network, errors="ignore").unlink()
        except FileNotFoundError:
            pass

    # Check if output_path attribute is available
    if not hasattr(path, "output_path"):
        raise AttributeError("Output path not set in config.yaml")
    sub_path = path.output_path

    # Check if output_path_private attribute is available
    #TODO: At the moment, need to have some private output path, even if not used
    if not hasattr(path, "output_path_private"):
        raise AttributeError("Private output path not set in config.yaml")
    sub_path_private = path.output_path_private
        
    out_pth_public = data_file_path("", network=network, sub_path=sub_path, errors="ignore")
    out_pth_private = data_file_path("", network=network, sub_path=sub_path_private, errors="ignore")

    if public:
        out_pth = out_pth_public
    else:
        out_pth = out_pth_private

    if delete:
        # Clear output directory, removing all files and subdirectories
        network, sub_path, files = data_file_list(network=network,
                                                  sub_path=sub_path,
                                                  errors="ignore")

        print(f'Deleting all files in {out_pth}')

        # If out_pth is a zip file, delete it
        if out_pth.suffix == ".zip" and out_pth.exists():
            out_pth.unlink()
        else:
            # For safety that out_pth is in data/network directory
            if out_pth.parents[1] == path.data and out_pth.parents[0].name == network:
                pass
            else:
                raise ValueError(f"{out_pth} must be in a data/network directory")

            # Delete all files in output directory
            for f in files:
                pth = out_pth / f
                if pth.is_file():
                    pth.unlink()
                elif pth.is_dir():
                    rmtree(pth)

    # If either out_pth is a zip file that doesn't, create them
    if out_pth_public.suffix == ".zip" and not out_pth_public.exists():
        with ZipFile(out_pth_public, "w") as f:
            pass
    if out_pth_private.suffix == ".zip" and not out_pth_private.exists():
        with ZipFile(out_pth_private, "w") as f:
            pass
    
    # Must run combined instruments first
    if combined:
        run_combined_instruments(network,
                                baseline=baseline, verbose=True,
                                monthly=monthly, species=species,
                                public=public)

    # If include is empty, process all instruments in release schedule
    if len(instrument_include) == 0:
        with open_data_file("data_release_schedule.xlsx", network=network) as frs:
            instruments = pd.ExcelFile(frs).sheet_names
    else:
        instruments = instrument_include

    # Processing
    for instrument in instruments:
        if instrument not in instrument_exclude:
            baseline_flag = {True: "git_pollution_flag", False: ""}[baseline]
            run_individual_instrument(network, instrument, 
                                    baseline=baseline_flag, verbose=True,
                                    monthly=monthly, species=species,
                                    public=public)

    # Incorporate README file into output directory or zip file
    try:
        readme_file = data_file_path(filename='README.md',
                                    network=network)
        copy_to_archive(readme_file, out_pth)
    except FileNotFoundError:
        print("No README file found")

    # If error log files have been created, warn the user
    if data_file_path("error_log_combined.txt", network=network, errors="ignore").exists():
        print("!!! Errors occurred during processing. See error_log_combined.txt for details")
    if data_file_path("error_log_individual.txt", network=network, errors="ignore").exists():
        print("!!! Errors occurred during processing. See error_log_individual.txt for details")


if __name__ == "__main__":

    print("####################################")
    print("#####Processing public archive######")
    print("####################################")
    run_all("agage", species = ["cf4"])

    # print("####################################")
    # print("#####Processing private archive#####")
    # print("####################################")
    # run_all("agage", species = ["ch4"], public=False)