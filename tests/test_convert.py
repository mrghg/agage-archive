import xarray as xr
import pandas as pd
import numpy as np

from agage_archive.convert import resample

def test_resample():
    """Test resample function"""

    # Create test dataset
    time = pd.date_range(start="2021-01-01", end="2021-01-02", freq="1min")
    data = np.random.rand(len(time))
    ds = xr.Dataset({"time": time, 
                    "mf": ("time", data),
                    "mf_repeatability": ("time", data * 0.1),
                    "sampling_period": ("time", np.ones(len(time)) * 60),
                    "baseline": ("time", np.ones(len(time)))})

    ds["mf"].attrs["units"] = "1e-9"
    ds["mf"].attrs["calibration_scale"] = "TU-87"

    ds["mf_repeatability"].attrs["units"] = "1e-9"
    ds["mf_repeatability"].attrs["calibration_scale"] = "TU-87"

    # Set the baseline variable to 0 once every two hours
    ds.baseline.values[::120] = 0

    # Test resample function
    ds_resample = resample(ds, resample_period=3600, resample_threshold=600)

    # Check that the dataset has been resampled to hourly
    assert ds_resample.time.diff("time").median() == np.timedelta64(3600, "s")

    # Check that the baseline variable has been resampled correctly
    assert np.all(ds_resample.baseline.values[::2] == 0)
    assert np.all(ds_resample.baseline.values[1::2] == 1)

    # Check that the sampling_period variable has been resampled correctly
    assert np.all(ds_resample.sampling_period.values == 3600)

    # Check that mf is the mean of the original data
    assert np.isclose(ds_resample.mf.values[0], data[:60].mean())

    # Check variability has been calculated correctly
    assert np.isclose(ds_resample.mf_variability.values[0], np.std(data[:60], ddof=1))

    # Check that units are consistent
    assert ds_resample.mf_variability.attrs["units"] == ds_resample.mf.attrs["units"]
    assert ds_resample.mf_variability.attrs["calibration_scale"] == ds_resample.mf.attrs["calibration_scale"]
    assert ds_resample.mf_repeatability.attrs["units"] == ds_resample.mf.attrs["units"]
