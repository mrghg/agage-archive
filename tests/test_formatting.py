import numpy as np
import xarray as xr
import pandas as pd
import warnings
from tempfile import NamedTemporaryFile
import json

from agage_archive.config import open_data_file, data_file_path
from agage_archive.run import run_individual_instrument


def check_cf_compliance(dataset):
    """Tests the CF compliance of a dataset when written to netCDF format.
    Taken from openghg/tests/helpers/cfchecking.py module

    Args: 
        dataset: An xarray.Dataset
    Returns:
        bool: True if compliant
    """
    from cfchecker.cfchecks import CFChecker

    checker = CFChecker(debug=False, version="1.8")

    catch_errors_warnings = False

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with NamedTemporaryFile(suffix=".nc") as tmpfile:
            dataset.to_netcdf(tmpfile.name)
            result = checker.checker(file=tmpfile.name)

            if result["global"]["FATAL"] or result["global"]["ERROR"]:
                catch_errors_warnings = True

            for var in result["variables"].keys():            
                if result["variables"][var]["FATAL"] or result["variables"][var]["ERROR"]:
                    catch_errors_warnings = True

    if catch_errors_warnings:
        return False
    else:
        return True


def test_cf_compliance():
    """Test CF compliance of the NF3 MHD dataset, which is generated by the run_single_site function on the agage-test data
    
    Args:
        None
    
    """
    
    network = 'agage_test'
    sub_path = 'output'
    instrument = "GCMS-Medusa"
    site = "MHD" # data_release_schedule modified so that this is the only site
    species = "nf3"

    # Get current version number from attributes.json
    with open_data_file("attributes.json", network=network) as f:
        attributes = json.load(f)
    version = attributes["version"]

    pth = data_file_path("",
                         network=network,
                         sub_path=sub_path,
                         errors="ignore")

    if not pth.exists():
        pth.mkdir()

    # Delete any files in pth
    for f in pth.rglob("*"):
        if f.is_file():
            f.unlink()

    run_individual_instrument(network=network,
                    instrument=instrument,
                    species=[species],
                    monthly=False,
                    verbose=False)
    
    with open_data_file(filename=f'{network.lower()}-{instrument.lower()}_{site}_{species}_{version}.nc',
                        network=network,
                        sub_path=sub_path + '/nf3/individual-instruments',
                        ) as f:
        ds = xr.load_dataset(f)

    assert check_cf_compliance(dataset=ds) == True
