import numpy as np
import xarray as xr
import pandas as pd
from agage_archive.formatting import monthly_baseline

def test_monthly_baseline():
    
    # Create a sample dataset
    time = pd.date_range("1991-01-01", "1991-12-31", freq="D")
    data = {"mf": (["time"], np.random.rand(len(time))),
            "mf_repeatability": (["time"], np.random.rand(len(time)))}
    ds = xr.Dataset(coords={"time": time}, data_vars=data)
    ds.mf.attrs["units"] = "1e-12"

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