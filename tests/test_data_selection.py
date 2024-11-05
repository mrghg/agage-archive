import numpy as np
import pandas as pd
import xarray as xr

from agage_archive.data_selection import read_data_exclude, calibration_scale_default


def test_calibration_scale_defaults():
    '''Test calibration_scale_default function'''

    assert calibration_scale_default("agage_test", "CO2") == "WMO-X2019"
    assert calibration_scale_default("agage_test", "CH4") == "TU-87"
    assert calibration_scale_default("agage_test", "CO2", scale_defaults_file="defaults-test") == "TESTING"


def test_read_data_exclude():
    '''Test read_data_exclude function'''

    # ALE CGO CH3CCl3 is flagged at 1978-10-12 06:00 - 07:00
    # Let's test that that point is removed, surrounding points are not
    ds = xr.Dataset(
        {
            "mf": (["time"], np.array([1., 2., 3.]))
        },
        coords = {
            "time": pd.date_range(start="1978-10-12 05:30", 
                                  end = "1978-10-12 07:30",
                                  periods=3)
        }
    )

    ds.attrs["network"] = "agage_test"

    ds = read_data_exclude(ds, "ch3ccl3", "CGO", "ALE")

    #This one should be nan
    assert np.isnan(ds["mf"][1].values)
    #These should not be nan
    assert not np.isnan(ds["mf"][0].values)
    assert ds["mf"][0].values == 1
    assert not np.isnan(ds["mf"][2].values)
    assert ds["mf"][2].values == 3

    # Test that the combined_only column works as expected
    # ALE CGO CH3CCl3 is flagged for combined only at 1972-01-01 06:00 - 07:00
    # Let's test that that point is removed, surrounding points are not, if a combined file is being created
    for combined in [True, False]:
        ds = xr.Dataset(
            {
                "mf": (["time"], np.array([1., 2., 3.]))
            },
            coords = {
                "time": pd.date_range(start="1972-01-01 05:30", 
                                    end = "1972-01-01 07:30",
                                    periods=3)
            }
        )

        ds.attrs["network"] = "agage_test"

        ds = read_data_exclude(ds, "ch3ccl3", "CGO", "ALE", combined=combined)

        #This one should be nan
        if combined:
            assert np.isnan(ds["mf"][1].values)
        else:
            assert not np.isnan(ds["mf"][1].values)
            assert ds["mf"][1].values == 2

        #These should not be nan
        assert not np.isnan(ds["mf"][0].values)
        assert ds["mf"][0].values == 1
        assert not np.isnan(ds["mf"][2].values)
        assert ds["mf"][2].values == 3

