import pandas as pd
import xarray as xr
import numpy as np

from agage_archive import Paths
from agage_archive.convert import scale_convert
from agage_archive.io import read_ale_gage, combine_datasets, open_data_file


paths = Paths(test=True)

# Scale conversion factors are used in multiple tests
scale_conversion = pd.read_csv(paths.root / "data/scale_convert.csv",
                        index_col="Species")


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

    ds_ale = read_ale_gage("cfc-11", "CGO", "ALE", testing_path=True)

    # test UTC conversion
    ds_ale_local = read_ale_gage("cfc-11", "CGO", "ALE", testing_path=True,
                                 utc=False, data_exclude=False)
    
    # check that all df_ale_local timestamps are 10 hours ahead of df_ale (CGO is UTC+10)
    assert ((ds_ale_local.time.to_series().index - ds_ale.time.to_series().index) == pd.Timedelta("10H")).all()

    # Check that some data have been excluded
    assert ds_ale.mf.to_series().isnull().sum() > ds_ale_local.mf.to_series().isnull().sum()

    ds_ale_noscale = read_ale_gage("cfc-11", "CGO", "ALE", testing_path=True,
                                scale=None)

    # Check that the scale conversion has been applied
    assert np.isclose(np.nanmean(ds_ale.mf.values / ds_ale_noscale.mf.values),
        scale_conversion.loc["cfc-11", "SIO-98/SIO-93"] * scale_conversion.loc["cfc-11", "SIO-05/SIO-98"])


def test_combine_datasets():

    species = "ch3ccl3"

    ds = combine_datasets(species, "CGO", testing_path=True)

    ds_ale_noscale = read_ale_gage(species, "CGO", "ALE", testing_path=True, scale=None)

    # Check that the scale conversion has been applied
    assert np.isclose(np.nanmean(ds.reindex_like(ds_ale_noscale).mf.values / ds_ale_noscale.mf.values),
        scale_conversion.loc[species, "SIO-98/SIO-93"] * scale_conversion.loc[species, "SIO-05/SIO-98"], rtol=0.00001)

    ds_gage_noscale = read_ale_gage(species, "CGO", "GAGE", testing_path=True, scale=None)

    # Check that the scale conversion has been applied
    assert np.isclose(np.nanmean(ds.reindex_like(ds_gage_noscale).mf.values / ds_gage_noscale.mf.values),
        scale_conversion.loc[species, "SIO-98/SIO-93"] * scale_conversion.loc[species, "SIO-05/SIO-98"], rtol=0.00001)


def test_open_data_file():

    #TODO: FINISH THIS

    # Open a global data file
    file_contents = open_data_file("attributes.json")

    # Open a network specific file
    file_contents = open_data_file("data-gcms-nc/AGAGE-GCMS-Medusa_CGO_ch3ccl3.nc",
                                   network = "agage_test")

    # Open a zipped file
    file_contents = open_data_file("...")