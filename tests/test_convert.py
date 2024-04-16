import xarray as xr
import pandas as pd
import numpy as np

from agage_archive.convert import resample
from agage_archive.convert import monthly_baseline


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
    ds["mf_repeatability"].attrs["long_name"] = "Repeatability"

    # Set the baseline variable to 0 once every two hours
    ds.baseline.values[::120] = 0

    # Test resample function
    ds_resample = resample(ds, resample_period="3600s", resample_threshold="600s")

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
    assert np.isclose(ds_resample.mf_variability.values[0], np.std(data[:60], ddof=0))

    # Check that units are consistent
    assert ds_resample.mf_variability.attrs["units"] == ds_resample.mf.attrs["units"]
    assert ds_resample.mf_variability.attrs["calibration_scale"] == ds_resample.mf.attrs["calibration_scale"]
    assert ds_resample.mf_repeatability.attrs["units"] == ds_resample.mf.attrs["units"]


def test_monthly_baseline():
    
    # Create a sample dataset
    time = pd.date_range("1991-01-01", "1991-12-31", freq="D")
    data = {"mf": (["time"], np.random.rand(len(time))),
            "mf_repeatability": (["time"], np.random.rand(len(time)))}
    ds = xr.Dataset(coords={"time": time}, data_vars=data)
    ds.mf.attrs["units"] = "1e-12"
    ds.mf.attrs["calibration_scale"] = "TU-87"
    ds.mf_repeatability.attrs["units"] = "1e-12"
    ds.mf_repeatability.attrs["calibration_scale"] = "TU-87"
    ds.mf_repeatability.attrs["long_name"] = "Repeatability"
    ds.attrs["version"] = "test"

    # Create a sample baseline dataset
    baseline_data = {"baseline": (["time"], np.random.choice([0, 1], size=len(time)))}
    ds_baseline = xr.Dataset(coords={"time": time}, data_vars=baseline_data)
    ds_baseline.attrs["baseline_flag"] = "Baseline Flag"

    ds_baseline_points = ds.where(ds_baseline.baseline == 1, drop=True)

    # Call the monthly_baseline function
    ds_monthly = monthly_baseline(ds, ds_baseline)

    # Check if the monthly mean is calculated correctly
    assert np.allclose(ds_monthly.mf.values, ds_baseline_points.mf.resample(time="1MS").mean().values)

    # Check if the monthly standard deviation is calculated correctly
    expected_std = ds_baseline_points.mf.resample(time="1MS").std()
    assert np.allclose(ds_monthly.mf_variability.values, expected_std.values)
    assert ds_monthly.mf_variability.attrs["long_name"] == "Monthly standard deviation of baseline mole fractions"
    assert ds_monthly.mf_variability.attrs["units"] == ds.mf.attrs["units"]

    # Check if the monthly standard error in mean is calculated correctly
    expected_se = ds_baseline_points.mf_repeatability.resample(time="1MS").mean() / np.sqrt(ds_baseline_points.mf.resample(time="1MS").count())
    assert np.allclose(ds_monthly.mf_repeatability, expected_se)
    assert ds_monthly.mf_repeatability.attrs["long_name"] == "Monthly standard error in mean of baseline mole fractions"
    assert ds_monthly.mf_repeatability.attrs["units"] == ds.mf.attrs["units"]

    # Check if the baseline flag is added correctly
    assert ds_monthly.attrs["baseline_flag"] == ds_baseline.attrs["baseline_flag"]

    # Check that a version number is present
    assert "version" in ds_monthly.attrs.keys()
