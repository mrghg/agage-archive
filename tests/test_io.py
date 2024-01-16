import pandas as pd
import xarray as xr
import numpy as np
import json

from agage_archive import Paths, open_data_file, data_file_path, data_file_list
from agage_archive.convert import scale_convert
from agage_archive.io import read_ale_gage, combine_datasets, read_nc_path, read_nc, read_baseline


paths = Paths("agage_test")


# Scale conversion factors are used in multiple tests
with open_data_file("scale_convert.csv") as f:
    scale_conversion = pd.read_csv(f, index_col="Species")


def test_scale_convert():

    def test_dataset(time, species, scale):
        coords = {"time": pd.date_range(time, time)}
        data = {"mf": (["time"], [1.]),
                "mf_repeatability": (["time"], [0.1])}
        ds = xr.Dataset(coords=coords, data_vars=data)
        ds.attrs["species"] = species
        ds.attrs["calibration_scale"] = scale
        return ds
    
    # Create xarray dataset with one time point
    ds = test_dataset("1991-01-01", "cfc-11", "SIO-93")

    # Test a set of conversion factors
    tests = [("cfc-11", "SIO-93", "SIO-98", "SIO-98/SIO-93"),
             ("cfc-11", "SIO-93", "SIO-05", "SIO-05/SIO-98*SIO-98/SIO-93"),
             ("ccl4", "SIO-93", "SIO-05", "SIO-05/SIO-98*SIO-98/SIO-93"),
             ("ch4", "CSIRO-94", "TU-87", "TU-87/CSIRO-94"),
             ("n2o", "SIO-93", "SIO-05", "SIO-05/SIO-98*SIO-98/SIO-93")]

    for test in tests:
        ds = test_dataset("1991-01-01", test[0], test[1])
        ds_new = scale_convert(ds, test[2])

        ratio = test[3].split("*")
        if len(ratio) == 1:
            assert ds_new.mf.values[0] / ds.mf.values[0] == \
                scale_conversion.loc[test[0], test[3]]
        else:
            assert ds_new.mf.values[0] / ds.mf.values[0] == \
                scale_conversion.loc[test[0], ratio[0]] * \
                scale_conversion.loc[test[0], ratio[1]]
                
    # Period where N2O time conversion applies
    ds = test_dataset("1985-01-01", "n2o", "SIO-93")
    ds_new = scale_convert(ds, "SIO-98")

    assert ds_new.mf.values[0] / ds.mf.values[0] == \
        scale_conversion.loc["n2o", "SIO-98/SIO-93"] * \
        0.9962230167482587

    # Period where N2O time conversion applies (inverse of above)
    ds = test_dataset("1985-01-01", "n2o", "SIO-98")
    ds_new = scale_convert(ds, "SIO-93")

    assert ds_new.mf.values[0] / ds.mf.values[0] == \
        1./(scale_conversion.loc["n2o", "SIO-98/SIO-93"] * \
        0.9962230167482587)

    # Test one where we go backwards
    ds = test_dataset("1991-01-01", "cfc-11", "SIO-05")
    ds_new = scale_convert(ds, "SIO-93")
    assert ds_new.mf.values[0] / ds.mf.values[0] == \
        (1./scale_conversion.loc["cfc-11", "SIO-98/SIO-93"]) * \
        (1./scale_conversion.loc["cfc-11", "SIO-05/SIO-98"])
    

def test_read_ale_gage():

    species = "ch3ccl3"

    ds_ale = read_ale_gage("agage_test", species, "CGO", "ALE")

    # test UTC conversion
    ds_ale_local = read_ale_gage("agage_test", species, "CGO", "ALE",
                                 utc=False, data_exclude=False)
    
    # check that all df_ale_local timestamps are 10 hours ahead of df_ale (CGO is UTC+10)
    assert ((ds_ale_local.time.to_series().index - ds_ale.time.to_series().index) == pd.Timedelta("10H")).all()

    # Check that some data have been excluded
    assert ds_ale.mf.to_series().isnull().sum() > ds_ale_local.mf.to_series().isnull().sum()

    ds_ale_noscale = read_ale_gage("agage_test", species, "CGO", "ALE",
                                scale=None)

    # Check that the scale conversion has been applied
    assert np.isclose(np.nanmean(ds_ale.mf.values / ds_ale_noscale.mf.values),
        scale_conversion.loc[species, "SIO-98/SIO-93"] * scale_conversion.loc[species, "SIO-05/SIO-98"])


def test_combine_datasets():

    network = "agage_test"

    species = "ch3ccl3"

    ds = combine_datasets(network, species, "CGO")

    ds_ale_noscale = read_ale_gage(network, species, "CGO", "ALE", scale=None)

    # Check that the scale conversion has been applied
    assert np.isclose(np.nanmean(ds.reindex_like(ds_ale_noscale).mf.values / ds_ale_noscale.mf.values),
        scale_conversion.loc[species, "SIO-98/SIO-93"] * scale_conversion.loc[species, "SIO-05/SIO-98"], rtol=0.00001)

    ds_gage_noscale = read_ale_gage(network, species, "CGO", "GAGE", scale=None)

    # Check that the scale conversion has been applied
    assert np.isclose(np.nanmean(ds.reindex_like(ds_gage_noscale).mf.values / ds_gage_noscale.mf.values),
        scale_conversion.loc[species, "SIO-98/SIO-93"] * scale_conversion.loc[species, "SIO-05/SIO-98"], rtol=0.00001)


def test_data_file_path():

    assert data_file_path("test.txt", "agage_test", "path_test_files").exists()

    # This should just return the zip file path, but checks for the existance of the file internally
    assert data_file_path("test_top_level.txt", "agage_test", "path_test_files/A.zip").exists()

    # This should just return the zip file path, but checks for the existance of the file internally
    assert data_file_path("B/C.txt", "agage_test", "path_test_files/A.zip").exists()

    # Test that we can find a required file in this repo
    assert data_file_path("attributes.json", this_repo=True).exists()


def test_open_data_file():

    # Test that the attributes json is openned correctly
    with open_data_file("attributes.json", this_repo=True) as f:
        attributes = json.load(f)
    assert "calibration_scale" in attributes.keys()
    assert attributes["species"] == ""

    # Test that we can open a file within the A.zip file
    with open_data_file("B/C.txt", "agage_test", "path_test_files/A.zip") as f:
        assert f.read().decode("utf-8") == "test"


def test_data_file_list():

    # Test that we can list files in a folder
    network, sub_path, files = data_file_list("agage_test", "path_test_files")
    assert network == "agage_test"
    assert sub_path == "path_test_files/"
    assert "test.txt" in files

    # Test that we can list files in a zip archive
    files = data_file_list("agage_test", "path_test_files/A.zip")[2]
    assert "test_top_level.txt" in files
    assert "B/C.txt" in files

    # Test that we can list files in a zip archive with pattern
    files = data_file_list("agage_test", "path_test_files/A.zip", pattern="*.txt")[2]
    assert "test_top_level.txt" in files
    assert "B/C.txt" in files

    # Test that we can list files within a subdirectory of a zip archive
    files_zip = data_file_list("agage_test", "path_test_files/A.zip", pattern="B/*.txt")[2]
    assert "test_top_level.txt" not in files_zip
    assert "B/C.txt" in files_zip

    # Test that we get the same output but from an unzipped directory
    files = data_file_list("agage_test", "path_test_files/A", "B/*.txt")[2]
    assert "test_top_level.txt" not in files
    assert "B/C.txt" in files


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

