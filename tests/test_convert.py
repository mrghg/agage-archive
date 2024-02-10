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

    # Test resample function
    ds_resample = resample(ds, resample_period=3600, resample_threshold=600)

    # Check that the dataset has been resampled to hourly
    assert ds_resample.time.diff("time").median() == np.timedelta64(3600, "s")

    # Check that the resampled dataset has the same number of time points as the original
    assert len(ds_resample.time) == len(ds.time) / 60

    # Check that the baseline variable has been resampled correctly
    #TODO: Improve this test
    assert np.all(ds_resample.baseline.values == 1)

    # Check that the sampling_period variable has been resampled correctly
    assert np.all(ds_resample.sampling_period.values == 3600)
