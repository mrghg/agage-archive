import pandas as pd
from shutil import rmtree
from zipfile import ZipFile
from multiprocessing import Pool
import time
import traceback
import os
import xarray as xr

from agage_archive.config import Paths, open_data_file, data_file_list, data_file_path, \
    copy_to_archive, output_path
from agage_archive.data_selection import read_release_schedule, read_data_combination
from agage_archive.io import combine_datasets, combine_baseline, \
    read_nc, read_baseline, read_ale_gage, read_gcwerks_flask, \
    output_dataset
from agage_archive.formatting import format_species
from agage_archive.convert import monthly_baseline
from agage_archive.definitions import instrument_number


# Set number of threads for multiprocessing
nthreads = 1


def get_error(e):
    """Get error message from exception

    Args:
        e (Exception): Exception object

    Returns:
        str: Error message
    """
    tb = traceback.extract_tb(e.__traceback__)

    stack_files_and_lines = []
    
    for t in tb:
        if "_archive" in t.filename:
            # Only include the filename and line no, not the full path
            stack_files_and_lines.append(f"{t.filename.split('/')[-1].split('.')[0]} (line {t.lineno})")

    error_type = type(e).__name__
    return f"{error_type} in stack: {' / '.join(stack_files_and_lines)}. {str(e)}"


def delete_archive(network, public = True):
    """Delete all files in output directory before running

    Args:
        network (str): Network for output filenames
        public (bool): Whether the dataset is for public release
    """

    path = Paths(network, errors="ignore", public=public)

    try:
        out_pth, _ = output_path(network, "_", "_", "_",
                                public=public,
                                errors="ignore_inputs")
    except FileNotFoundError:
        print(f"Output directory or archive for does not exist, continuing")
        return

    # Find all files in output directory
    _, _, files = data_file_list(network=network,
                                sub_path=path.output_path,
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


def create_empty_archive(network, public = True):
    """Create an empty output zip file or folders

    Args:
        network (str): Network for output filenames
        public (bool): Whether the dataset is for public release
    """

    out_pth, _ = output_path(network, "_", "_", "_",
                            public=public,
                            errors="ignore")

    if out_pth.suffix == ".zip" and not out_pth.exists():
        with ZipFile(out_pth, "w") as f:
            pass
    elif out_pth.is_dir():
        out_pth.mkdir(parents=True, exist_ok=True)


def run_timestamp_checks(ds,
                        ds_baseline=None,
                        species="",
                        site=""):

    # Check for duplicate time stamps
    timestamps = ds["time"].to_series()
    if timestamps.duplicated().any():
        # Create list of duplicated timestamps
        duplicated = timestamps[timestamps.duplicated()].unique()
        duplicated_str = ", ".join([str(d) for d in duplicated])

        # List of instrument types that have duplicate timestamps
        instrument_types = ds["instrument_type"].to_series()
        instrument_types = instrument_types[timestamps.duplicated()].unique()

        # find instrument name in instrument_number
        instrument_names = [k for k, v in instrument_number.items() if v in instrument_types]
        instrument_names = ", ".join(instrument_names)

        raise ValueError(f"Duplicate timestamps in {species} at {site}: {duplicated_str} for instrument {instrument_names}")

    if ds_baseline:
        if ds_baseline["time"].to_series().duplicated().any():
            raise ValueError(f"Duplicate timestamps in baseline for {species} at {site}")
    
    # check that the time stamps are the same in the data and baseline files
    if ds_baseline:
        if (ds_baseline.time != ds.time).any():
            raise ValueError(f"Data and baseline files for {species} at {site} have different timestamps")
        

def run_individual_site(site, species, network, instrument,
                        rs, read_function, read_baseline_function, instrument_out,
                        baseline=False,
                        monthly=False,
                        verbose=False,
                        public=True,
                        resample=True,
                        top_level_only=False):
    """Process individual data files for a given site.
    Reads the release schedule for the site

    Args:
        site (str): Site to process. Must match sheet names in release schedule, e.g.:
            "MHD", "MLO", "SMO", ...
        species (str): species to process. If empty, process all species
        network (str): Network for output filenames
        instrument (str): Instrument to process. Must match sheet names in release schedule, e.g.:
            "AGAGE", "ALE", "GAGE", "GCMD", ...
        rs (pd.DataFrame): Release schedule
        read_function (function): Function to read data files
        read_baseline_function (function): Function to read baseline files
        instrument_out (str): Instrument name for output filenames
        baseline (bool): Process baselines. Boolean as only one baseline flag is available (GIT)
        monthly (bool): Produce monthly baseline files
        verbose (bool): Print progress to screen
        public (bool, optional): Whether the dataset is for public release. Default to True.
        resample (bool, optional): Whether to resample the data, if needed. Default to True.
        top_level_only (bool, optional): Whether to only output to the top-level directory, 
            and ignore the individual instrument folder. Default to False.
    """

    if read_function.__name__ == "read_gcwerks_flask":
        site_str = site.lower()
    else:
        site_str = ""

    paths = Paths(network, public=public, errors="ignore_outputs", site = site_str)

    error_log = []

    try:

        if rs.loc[species, site].lower() != "x":

            ds = read_function(network, species, site, instrument,
                            public = public, verbose=verbose, resample=resample)

            if baseline:
                ds_baseline = read_baseline_function(network, species, site, instrument,
                                            flag_name = baseline,
                                            verbose = verbose)
            else:
                ds_baseline = None
                
            run_timestamp_checks(ds, ds_baseline, species, site)

            # If multiple instruments, store individual file in subdirectory
            instrument_dates = read_data_combination(network, species, site,
                                                    verbose=False)

            if top_level_only:
                folders = []
            else:
                folders = [f"{species}/individual-instruments"]

            # if there is no combined data file, also store individual file in top-level directory
            if len(instrument_dates) <= 1:
                # Add to top-level directory
                folders.append(f"{species}")
            else:
                if top_level_only:
                    raise ValueError(f"Looks like combined instruments has been run for {species} at {site}, but top_level_only is set to True")
            

            for output_subpath in folders:

                if "individual" in output_subpath:
                    instrument_str = instrument_out
                else:
                    instrument_str = ""

                # Check if file already exists in the top-level directory. 
                # This can happen if only one instrument is specified in data_combination.xlsx
                if output_subpath == f"{species}":
                    if data_file_list(network=network,
                                    sub_path=paths.output_path,
                                    pattern = f"{format_species(species)}/{network.lower()}_{site.lower()}_{format_species(species)}*.nc",
                                    errors="ignore")[2]:
                        raise FileExistsError(f"Top-level file already exists for {species} at {site} (now trying to add instrument {instrument_out}). "\
                                            "Add to data_combination.xlsx to tell me how to combine the data.")

                output_dataset(ds, network, instrument=instrument_str,
                            output_subpath=output_subpath,
                            end_date=rs.loc[species, site],
                            public=public,
                            verbose=verbose)

                if baseline:
                    if (ds_baseline.time != ds.time).any():
                        raise ValueError(f"Baseline and data files for {species} at {site} have different timestamps")
                    output_dataset(ds_baseline, network, instrument=instrument_str,
                            output_subpath=output_subpath + "/baseline-flags",
                            end_date=rs.loc[species, site],
                            extra="git-baseline",
                            public=public,
                            verbose=verbose)
                    
                    if monthly:
                        ds_baseline_monthly = monthly_baseline(ds, ds_baseline)
                        output_dataset(ds_baseline_monthly, network, instrument=instrument_str,
                            output_subpath=output_subpath + "/monthly-baseline",
                            end_date=rs.loc[species, site],
                            extra="monthly-baseline",
                            public=public,
                            verbose=verbose)

                else:
                    if monthly:
                        raise NotImplementedError("Monthly baseline files can only be produced if baseline flag is specified")

            error_log.append("")

        else:

            error_log.append("")

    except Exception as e:

        error_log.append(get_error(e))
    
    return (site, species, error_log[0])


def run_individual_instrument(network, instrument,
                              verbose = False,
                              baseline = "",
                              monthly = False,
                              species = [],
                              public=True,
                              resample=True,
                              top_level_only=False):
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
        resample (bool, optional): Whether to resample the data, if needed. Default to True.
        top_level_only (bool, optional): Whether to only output to the top-level directory,
            and ignore the individual instrument folder. Default to False.
    """
    
    rs = read_release_schedule(network, instrument, public=public)

    if instrument.upper() == "ALE" or instrument.upper() == "GAGE":
        read_function = read_ale_gage
        read_baseline_function = read_baseline
        instrument_out = instrument.lower() + "-gcmd"
    elif instrument.upper() == "GCMS-MEDUSA-FLASK":
        read_function = read_gcwerks_flask
        read_baseline_function = None
        instrument_out = "gcms-medusa-flask"
    else:
        read_function = read_nc
        read_baseline_function = read_baseline
        instrument_out = instrument.lower()

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

        with Pool(processes=nthreads) as p:
            result = p.starmap_async(run_individual_site, [(site, sp, network, instrument,
                                                            rs, read_function, read_baseline_function, instrument_out,
                                                            baseline, monthly, verbose, public, resample, top_level_only) for site in rs.columns])

            for r in result.get():
                error_log.append(r)

    has_errors = any([error[2] for error in error_log])

    if has_errors:
        # save errors to file
        with open(data_file_path("error_log_individual.txt", network=network, errors="ignore"), "w") as f:
            # write the date and time of the error
            f.write("Processing attempted on " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
            for error in error_log:
                if error[2]:
                    f.write(f"{error[0]} {error[1]}: {error[2]}\n")


def run_combined_site(site, species, network, 
                    baseline=False,
                    monthly=False,
                    verbose=False,
                    public=True,
                    resample=True):
    """Process combined data files for a given site.
    Reads the data selection file to determine which species to process

    Args:
        site (str): Site to process. Must match sheet names in data selection file
        species (list): List of species to process. If empty, process all species
        network (str): Network for output filenames
        baseline (bool): Process baselines. Boolean as only one baseline flag is available (GIT)
        monthly (bool): Produce monthly baseline files
        verbose (bool): Print progress to screen
        public (bool, optional): Whether the dataset is for public release. Default to True.
        resample (bool, optional): Whether to resample the data, if needed. Default to True.
    """

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
            return [(site, "None", "")]
    else:
        # Process all species in the data selection file
        species_to_process = df.index.values

    error_log = []

    # Loop through species in index
    for sp in species_to_process:

        try:

            # Produce combined dataset
            if verbose:
                print(f"... combining datasets for {sp} at {site}")
            ds = combine_datasets(network, sp, site,
                                verbose=verbose, public=public, resample=resample)

            if baseline:
                if verbose:
                    print(f"... combining baselines for {sp} at {site}")
                # Note that GIT baselines is hard-wired here because Met Office not available for ALE/GAGE
                ds_baseline = combine_baseline(network, sp, site,
                                            verbose=verbose, public=public)

            else:
                ds_baseline = None

            # Check for duplicate time stamps
            run_timestamp_checks(ds, ds_baseline, sp, site)

            output_subpath = f"{sp}"

            if verbose:
                print(f"... outputting combined dataset for {sp} at {site}")
            output_dataset(ds, network,
                        output_subpath=output_subpath,
                        instrument="",
                        public=public,
                        verbose=verbose)
            
            if baseline:
                if verbose:
                    print(f"... outputting combined baseline for {sp} at {site}")
                output_dataset(ds_baseline, network,
                            output_subpath=output_subpath + "/baseline-flags",
                            instrument="",
                            extra="git-baseline",
                            public=public,
                            verbose=verbose)

                if monthly:
                    ds_baseline_monthly = monthly_baseline(ds, ds_baseline)
                    output_dataset(ds_baseline_monthly, network,
                            output_subpath=output_subpath + "/monthly-baseline",
                            instrument="",
                            extra="monthly-baseline",
                            public=public,
                            verbose=verbose)

            else:
                if monthly:
                    raise NotImplementedError("Monthly baseline files can only be produced if baseline flag is specified")

            error_log.append("")

        except Exception as e:

            error_log.append(get_error(e))

    return [(site, sp, error) for sp, error in zip(species_to_process, error_log)]


def run_combined_instruments(network,
                             baseline = False,
                             monthly = False,
                             verbose = False,
                             species = [],
                             public = True,
                             resample=True):
    """Process combined data files for a given network.
    Reads the data selection file to determine which sites to process

    Args:
        network (str): Network for output filenames
        baseline (bool): Process baselines. Boolean as only one baseline flag is available (GIT)
        monthly (bool): Produce monthly baseline files
        verbose (bool): Print progress to screen
        species (list): List of species to process. If empty, process all species
        public (bool, optional): Whether the dataset is for public release. Default to True.
        resample (bool, optional): Whether to resample the data, if needed. Default to True.
    """

    if not isinstance(species, list):
        raise TypeError("Species must be a list")

    with open_data_file("data_combination.xlsx", network=network) as data_selection_path:
        sites = pd.ExcelFile(data_selection_path).sheet_names

    with Pool(processes=nthreads) as p:
        args = [(site, species, network, baseline, monthly, verbose, public, resample) for site in sites]
        result = p.starmap_async(run_combined_site, [(site, species, network, baseline, monthly, verbose, public, resample) for site in sites])
    
        error_log = []

        for r in result.get():
            error_log.extend(r)

    has_errors = any([error[2] for error in error_log])

    if has_errors:
        # save errors to file
        with open(data_file_path("error_log_combined.txt", network=network, errors="ignore"), "w") as f:
            # write the date and time of the error
            f.write("Processing attempted on " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
            for error in error_log:
                if error[2]:
                    f.write(f"{error[0]} {error[1]}: {error[2]}\n")


def run_all(network,
            delete = True,
            combined = True,
            baseline = True,
            monthly = True,
            instrument_include = [],
            instrument_exclude = ["GCPDD"],
            species = [],
            public = True,
            resample=True,
            top_level_only=False,):
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
        resample (bool, optional): Whether to resample the data, if needed. Default to True.
        top_level_only (bool, optional): Whether to only output to the top-level directory,
            and ignore the individual instrument folder. Default to False.
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

    path = Paths(network, public = public, errors="ignore")

    # Delete log files, if they exist
    for log_file in ["error_log_combined.txt", "error_log_individual.txt"]:
        try:
            data_file_path(log_file, network=network, errors="ignore").unlink()
        except FileNotFoundError:
            pass

    # Check if output_path attribute is available
    if not hasattr(path, "output_path"):
        raise AttributeError("Output path not set in config.yaml")

    if delete:
        delete_archive(network, public=public)
        
    # If either out_pth is a zip file that doesn't exist, create
    create_empty_archive(network, public=public)

    # Must run combined instruments first
    if combined:
        run_combined_instruments(network,
                                baseline=baseline, verbose=True,
                                monthly=monthly, species=species,
                                public=public, resample=resample)

    # If include is empty, process all instruments in release schedule
    if len(instrument_include) == 0:
        with open_data_file("data_release_schedule.xlsx", network=network, errors = "ignore") as frs:
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
                                    public=public, resample=resample, top_level_only=top_level_only)

    # Incorporate README file into output directory or zip file
    try:
        readme_file = data_file_path(filename='README.md',
                                    network=network, errors = "ignore_inputs")
        copy_to_archive(readme_file, network, public=public)
    except FileNotFoundError:
        print("No README file found")

    # If error log files have been created, warn the user
    if data_file_path("error_log_combined.txt", network=network, errors="ignore").exists():
        print("!!! Errors occurred during processing. See error_log_combined.txt for details")
    if data_file_path("error_log_individual.txt", network=network, errors="ignore").exists():
        print("!!! Errors occurred during processing. See error_log_individual.txt for details")


def read_repeatability(std_file, species):
    """Read the repeatability values from the standard file
    """

    # Read the second and third lines of the file
    with open(std_file) as f:
        header = f.readline()
        colnames1 = f.readline()
        colnames2 = f.readline()

    colnames = [(col1 + col2).replace("-", "") for col1, col2 in  zip(colnames1.split(), colnames2.split())]
    colnames = [col[:-1].lower() if col[-1] == "C" else col for col in colnames]

    df = pd.read_csv(std_file,
                skiprows=3, delim_whitespace=True,
                names=colnames,
                parse_dates={'datetime': ['date', 'time']},
                date_format='%y%m%d %H%M',
                na_values="nan",
                index_col='datetime')[species.lower()]

    # Take running standard deviation with 48 hour window. Ignore NaNs
    df = df.rolling(window=24, min_periods=1, center=True).std()

    return df


def replace_repeatability_and_height(raw_file, std_file, species, site):
    """Replace the repeatability values in the raw file with the new values from the std file
    Also replace the inlet height values for TAC and TOB

    """

    with xr.open_dataset(raw_file, engine="netcdf4") as f:
        ds = f.load()

    # Read new repeatability file
    repeatability_df = read_repeatability(std_file, species)
    original_repeatability = ds['mf_repeatability'].to_series()

    # Replace the values
    ds['mf_repeatability'] = repeatability_df.reindex(original_repeatability.index, method='nearest')

    if(site == "TAC"):
        ds['inlet_height'].loc[ds["time"] < pd.to_datetime("2017-03-09 15:00")] = 100
        ds['inlet_height'].loc[ds["time"] >= pd.to_datetime("2017-03-09 15:00")] = 185
    elif(site == "TOB"):
        ds["inlet_height"].values[:] = 12

    # Save the new file
    ds.to_netcdf(raw_file)


def preprocess():
    """Preprocess data files before running the main script
    Hopefully this will become redundant in the future
    
    """

    paths = Paths("agage")

    # SOME MD DATA IS MISSING AND IN THE WRONG FORMAT
    sites = {"TAC": ["n2o", "co"],
            "TOB": ["sf6"]}
    og_paths = {"TAC": "/agage/summary/netcdf-decc/md/AGAGE-GCMD_TAC_species.nc",
                "TOB": "/agage/summary/netcdf-other/md/AGAGE-GCECD_TOB_species.nc"}

    for site in sites:

        for sp in sites[site]:

            # Copy MD files from md to md-modified directory, and replace repeatability
            md_folder = data_file_path("", "agage", sub_path=paths.md_path)
            og_path = og_paths[site].replace('species', sp.lower())
            og_file = og_path.split('/')[-1]
            std_path = og_path.replace(og_file, "") + f"stds/{site.upper()}_GCMD_stds.dat"

            # Copy from to md path
            os.system(f"cp {og_path} {md_folder}")

            # Replace repeatability
            replace_repeatability_and_height(md_folder / og_file, std_path, sp, site)

    # For CGO H2 data, the PDD is mis-labelled (Issue #47)
    if (md_folder / "AGAGE-GCMD_CGO_h2_pdd.nc").exists():
        os.system(f"cp {md_folder / 'AGAGE-GCMD_CGO_h2_pdd.nc'} {md_folder / 'AGAGE-GCPDD_CGO_h2.nc'}")


if __name__ == "__main__":

    preprocess()

    start_time = time.time()

    print("####################################")
    print("#####Processing public archive######")
    print("####################################")
    run_all("agage", species = ["cfc-11"], public=True)

    # print("####################################")
    # print("#####Processing private archive#####")
    # print("####################################")
    # run_all("agage", species = ["cfc-11"], public=False)

    print(f"Time taken: {time.time() - start_time:.2f} seconds")