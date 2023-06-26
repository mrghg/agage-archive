from agage_archive import io
import pandas as pd
import xarray as xr

paths = io.Paths()

def test_scale_convert():

    def test_dataset(time, species, scale):
        coords = {"time": pd.date_range(time, time)}
        data = {"mf": (["time"], [1.])}
        ds = xr.Dataset(coords=coords, data_vars=data)
        ds.attrs["species"] = species
        ds.attrs["calibration_scale"] = scale
        return ds
    
    scale_conversion = pd.read_csv(paths.root / "data/scale_convert.csv",
                                   index_col="Species")

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
        ds_new = io.scale_convert(ds, test[2])

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
    ds_new = io.scale_convert(ds, "SIO-05")

    assert ds_new.mf.values[0] / ds.mf.values[0] == \
        scale_conversion.loc["n2o", "SIO-98/SIO-93"] * \
        0.9962230167482587

