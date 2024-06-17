import pandas as pd
import xarray as xr
import numpy as np

from agage_archive.config import Paths, open_data_file
from agage_archive.io import read_ale_gage, read_nc, combine_datasets, read_nc_path, \
    read_baseline, combine_baseline, output_dataset, read_gcwerks_flask
from agage_archive.convert import scale_convert


paths = Paths("agage_test")

def test_read_ale_gage():

    species = "ch3ccl3"

    ds_ale = read_ale_gage("agage_test", species, "CGO", "ALE")

    # test UTC conversion
    ds_ale_local = read_ale_gage("agage_test", species, "CGO", "ALE",
                                 utc=False, data_exclude=False)
    
    # check that all df_ale_local timestamps are 10 hours ahead of df_ale (CGO is UTC+10)
    assert ((ds_ale_local.time.to_series().index - ds_ale.time.to_series().index) == pd.Timedelta("10h")).all()

    # Check that some data have been excluded
    assert ds_ale.mf.to_series().isnull().sum() > ds_ale_local.mf.to_series().isnull().sum()

    ds_ale_noscale = read_ale_gage("agage_test", species, "CGO", "ALE",
                                scale=None)

    # Check that the scale conversion has been applied (ch3ccl3 factors from openghg_calscales)
    assert np.isclose(np.nanmean(ds_ale.mf.values / ds_ale_noscale.mf.values), 0.9957*1.0184)


def test_combine_datasets():

    network = "agage_test"

    species = "ch3ccl3"

    ds = combine_datasets(network, species, "CGO")

    ds_ale_noscale = read_ale_gage(network, species, "CGO", "ALE", scale=None)

    # Check that the scale conversion has been applied
    assert np.isclose(np.nanmean(ds.reindex_like(ds_ale_noscale).mf.values / ds_ale_noscale.mf.values),
        0.9957*1.0184, rtol=0.00001)

    ds_gage_noscale = read_ale_gage(network, species, "CGO", "GAGE", scale=None)

    # Check that the scale conversion has been applied
    assert np.isclose(np.nanmean(ds.reindex_like(ds_gage_noscale).mf.values / ds_gage_noscale.mf.values),
        0.9957*1.0184, rtol=0.00001)

    # Test that a version number has been added
    assert "version" in ds.attrs.keys()
    assert ds.attrs["version"] != ""


def test_read_nc_path():

    nc_file, sub_path = read_nc_path("agage_test", "CH3CCl3", "CGO", "GCMS-Medusa")
    assert nc_file == "AGAGE-GCMS-Medusa_CGO_ch3ccl3.nc"
    assert sub_path == "data-gcms-nc"

    nc_file, sub_path = read_nc_path("agage_test", "NF3", "MHD", "GCMS-Medusa")
    assert nc_file == "AGAGE-GCMS-Medusa_MHD_nf3.nc"
    assert sub_path == "data-gcms-nc"

    nc_file, sub_path = read_nc_path("agage_test", "CH3CCl3", "CGO", "GCMD")
    assert nc_file == "AGAGE-GCMD_CGO_ch3ccl3.nc"
    assert sub_path == "data-nc"


def test_read_baseline():

    for flag_name in ["git_pollution_flag"]:
#    for flag_name in ["met_office_baseline_flag", "git_pollution_flag"]:

        ds_baseline = read_baseline("agage_test",
            "CH3CCl3", "CGO", "GCMS-Medusa",
            flag_name = flag_name)

        # Check baseline flags exist and has integer values
        assert "baseline" in ds_baseline.data_vars.keys()
        assert isinstance(ds_baseline.baseline.values[0], np.int8)

        # Check some attributes
        assert ds_baseline.attrs["species"] == "ch3ccl3"
        assert ds_baseline.attrs["instrument"] == "GCMS-Medusa"
        assert ds_baseline.attrs["site_code"] == "CGO"
        assert "citation" in ds_baseline.attrs.keys()

        assert ds_baseline.time.attrs["long_name"] == "time"

        # Check that the right baseline flag has been extracted
        if flag_name == "met_office_baseline_flag":
            assert "Met Office" in ds_baseline.attrs["comment"]
        elif flag_name == "git_pollution_flag":
            assert "Georgia Tech" in ds_baseline.attrs["comment"]

        # Test outputting baseline dataset
        output_dataset(ds_baseline, "agage_test", instrument="GCMS-Medusa",
            output_subpath="baselines/",
            extra = "-git-baseline",
            verbose=False)
        

def test_combine_baseline():

    ds_baseline = combine_baseline("agage_test", "CH3CCl3", "CGO")

    # Check baseline flags exist and has integer values
    assert "baseline" in ds_baseline.data_vars.keys()
    assert isinstance(ds_baseline.baseline.values[0], np.int8)

    # Check some attributes
    assert ds_baseline.attrs["species"] == "ch3ccl3"
    assert ds_baseline.attrs["site_code"] == "CGO"
    assert "citation" in ds_baseline.attrs.keys()

    assert ds_baseline.time.attrs["long_name"] == "time"

    # Check that ds_baseline has the same timestamps as the output of combine_datasets
    ds = combine_datasets("agage_test", "CH3CCl3", "CGO", verbose=False)
    assert (ds_baseline.time.values == ds.time.values).all()

    # Test that a version number has been added
    assert "version" in ds_baseline.attrs.keys()
    assert ds_baseline.attrs["version"] != ""


def test_timestamp():
    """Test that some timestamps are read correctly and corrected for the sampling time offset
    """

    # Medusa files have a substantial sampling period
    filename, sub_path = read_nc_path("agage_test", "CH3CCl3", "CGO", "GCMS-Medusa")

    with open_data_file(filename, network="agage_test", sub_path=sub_path) as f:
        ds_original = xr.open_dataset(f).load()
        
    ds = read_nc("agage_test", "CH3CCl3", "CGO", "GCMS-Medusa")

    time_offset = int(ds_original.time.attrs["sampling_time_seconds"])
    assert np.all(ds_original.time.values == ds.time.values + pd.Timedelta(seconds=time_offset)/2)

    assert ds.time.attrs["long_name"] == "time"
    assert ds.time.attrs["comment"] == "Timestamp is the start of the sampling period in UTC"


def test_picarro():
    """Test that Picarro file is read correctly, resampled to hourly and that the timestamp is at the beginning of the sampling period
    """

    filename, sub_path = read_nc_path("agage_test", "ch4", "THD", "Picarro")
    
    # Original data
    with open_data_file(filename, network="agage_test", sub_path=sub_path) as f:
        ds_original = xr.open_dataset(f).load()
    
    # Output of read_nc
    ds = read_nc("agage_test", "ch4", "THD", "Picarro")

    # Apply any scale conversion
    ds_original = scale_convert(ds_original, ds.attrs["calibration_scale"])

    # Find the first hour of the dataset
    first_hour = ds_original.time.dt.floor("h").min()
    first_hour_mean = ds_original.mf.sel(time=slice(first_hour, first_hour + pd.Timedelta("1h"))).mean()

    # Check that the dataset has been resampled to hourly
    assert np.isclose(ds.mf.values[0], first_hour_mean.values, rtol=0.0001)

    # Check that timestamp is at the beginning of the sampling period
    assert ds.time[0].values == first_hour.values

    # Check that the time variable has the correct attributes
    assert ds.time.attrs["long_name"] == "time"
    assert "Timestamp is the start of the sampling period in UTC" in ds.time.attrs["comment"]
    assert "Resampled" in ds.time.attrs["comment"]  


def test_read_gcwerks_flask():

    ds = read_gcwerks_flask("agage_test", "cf4", "CBW", "GCMS-Medusa-flask")

    # Check that the dataset has the correct attributes
    assert ds.attrs["species"] == "cf4"

    # First sample should be at 2021-02-15 1300 minus 30 minutes
    assert ds.time[0].values == pd.Timestamp("2021-02-15 1230")

    # Check that the time variable has the correct attributes
    assert ds.time.attrs["long_name"] == "time"
    assert "Timestamp is the start of sampling period" in ds.time.attrs["comment"]

    # Check that some attributes that are only in site attributes have been added
    assert "sampling_period" in ds.attrs.keys()
    assert ds.attrs["inlet_latitude"] == 2

