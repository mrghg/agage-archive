import numpy as np
import xarray as xr
import pandas as pd
import warnings
from tempfile import NamedTemporaryFile
from agage_archive.formatting import monthly_baseline
from agage_archive import open_data_file
from agage_archive.data_selection import read_release_schedule, read_data_combination
from agage_archive.io import read_nc, read_baseline, read_ale_gage, output_dataset

def test_monthly_baseline():
    
    # Create a sample dataset
    time = pd.date_range("1991-01-01", "1991-12-31", freq="D")
    data = {"mf": (["time"], np.random.rand(len(time))),
            "mf_repeatability": (["time"], np.random.rand(len(time)))}
    ds = xr.Dataset(coords={"time": time}, data_vars=data)
    ds.mf.attrs["units"] = "1e-12"
    ds.attrs["version"] = "test"

    # Create a sample baseline dataset
    baseline_data = {"baseline": (["time"], np.random.choice([0, 1], size=len(time)))}
    ds_baseline = xr.Dataset(coords={"time": time}, data_vars=baseline_data)
    ds_baseline.attrs["baseline_flag"] = "Baseline Flag"

    ds_baseline_points = ds.where(ds_baseline.baseline == 1, drop=True)

    # Call the monthly_baseline function
    ds_monthly = monthly_baseline(ds, ds_baseline)

    # Check if the monthly mean is calculated correctly
    assert (ds_monthly.mf == ds_baseline_points.resample(time="1MS").mean().mf).all()

    # Check if the monthly standard deviation is calculated correctly
    expected_std = ds_baseline_points.mf.resample(time="1MS").std()
    assert (ds_monthly.mf_variability == expected_std).all()
    assert ds_monthly.mf_variability.attrs["long_name"] == "Monthly standard deviation of baseline mole fractions"
    assert ds_monthly.mf_variability.attrs["units"] == ds.mf.attrs["units"]

    # Check if the monthly standard error in mean is calculated correctly
    expected_se = ds_baseline_points.resample(time="1MS").mean().mf_repeatability / np.sqrt(ds_baseline_points.mf.resample(time="1MS").count())
    assert (ds_monthly.mf_repeatability == expected_se).all()
    assert ds_monthly.mf_repeatability.attrs["long_name"] == "Monthly standard error in mean of baseline mole fractions"
    assert ds_monthly.mf_repeatability.attrs["units"] == ds.mf.attrs["units"]

    # Check if the baseline flag is added correctly
    assert ds_monthly.attrs["baseline_flag"] == ds_baseline.attrs["baseline_flag"]

    # Check that a version number is present
    assert "version" in ds_monthly.attrs.keys(

def check_cf_compliance(dataset):
    """Tests the CF compliance of a dataset when written to netCDF format.
    Taken from openghg/tests/helpers/cfchecking.py module

    Args: 
        dataset: An xarray.Dataset
    Returns:
        bool: True if compliant
    """
    from cfchecker import CFChecker

    checker = CFChecker(debug=False, version="1.8")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with NamedTemporaryFile(suffix=".nc") as tmpfile:
            dataset.to_netcdf(tmpfile.name)
            result = checker.checker(file=tmpfile.name)

            fatal_global = result["global"]["FATAL"]
            error_global = result["global"]["ERROR"]

            fatal_vars = []
            error_vars = []
            for var in result["variables"].keys():
                fatal_vars.append(result["variables"][var]["FATAL"])
                error_vars.append(result["variables"][var]["ERROR"])
            
            if fatal_global or error_global or fatal_vars or error_vars:
                return False
            else:
                return True

def run_single_site(network, instrument, site,
                              verbose = False,
                              baseline = "",
                              monthly = False,
                              species = [],
                              ):
    """Process individual data files for a given instrument at a single site.
    Reads the release schedule for the instrument

    Args:
        instrument (str): Instrument to process. Must match sheet names in release schedule, e.g.:
            "AGAGE", "ALE", "GAGE", "GCMD", ...
        verbose (bool): Print progress to screen
        baseline (str): Baseline flag to use. If empty, don't process baselines
        monthly (bool): Produce monthly baseline files
        species (list): List of species to process. If empty, process all species
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

    if species:
        # Process only those species that are in the release schedule
        species_to_process = [sp for sp in species if sp in rs.index.values]
        if not species_to_process:
            print(f"No species to process for {instrument}, skipping...")
            return
    else:
        # Process all species in the release schedule
        species_to_process = rs.index.values

    # Process for all species and sites
    for sp in species_to_process:
        
        if rs.loc[sp, site].lower() != "x":

            ds = read_function(network, sp, site, instrument,
                            verbose=verbose)

            if baseline:
                ds_baseline = read_baseline_function(network, sp, site, instrument,
                                            flag_name = baseline,
                                            verbose = verbose)

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
                            verbose=verbose)

            if baseline:
                if (ds_baseline.time != ds.time).any():
                    raise ValueError(f"Baseline and data files for {sp} at {site} have different timestamps")
                output_dataset(ds_baseline, network, instrument=instrument_out,
                            output_subpath=output_subpath + "/baseline_flags",
                            end_date=rs.loc[sp, site],
                            extra="-git-baseline",
                            verbose=verbose)
                
                if monthly:
                    ds_baseline_monthly = monthly_baseline(ds, ds_baseline)
                    output_dataset(ds_baseline_monthly, network, instrument=instrument_out,
                            output_subpath=output_subpath.replace("event", "monthly"),
                            end_date=rs.loc[sp, site],
                            extra="-monthly",
                            verbose=verbose)

            else:
                if monthly:
                    raise NotImplementedError("Monthly baseline files can only be produced if baseline flag is specified")   

def test_cf_compliance_MHD_NF3():
    """Test CF compliance of the NF3 MHD dataset, which is generated by the run_single_site function on the agage-test data
    
    Args:
        None
    
    """
    
    run_single_site(network='agage_test',
                    instrument='GCMS-Medusa',
                    site='MHD',
                    species=['nf3]'])
    
    with open_data_file(filename='AGAGE_TEST-GCMS-MEDUSA_MHD_nf3.nc',
                        network='agage_test',
                        sub_path='output/event/nf3'
                        ) as f:
        ds = xr.load_dataset(f)
    assert check_cf_compliance(dataset=ds)
