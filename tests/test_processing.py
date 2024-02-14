from agage_archive.config import Paths, open_data_file
from agage_archive.formatting import lookup_locals_and_attrs, \
    format_species, format_units, format_calibration_scale
from agage_archive.data_selection import calibration_scale_default


path = Paths()

def test_lookup_locals_and_attrs():
    local = {"species": None}
    attrs = {"species": "CFC-11"}
    assert lookup_locals_and_attrs("species", local, attrs) == format_species("CFC-11")

    local = {"species": "CFC-12"}
    attrs = {"species": "BLAH"}
    assert lookup_locals_and_attrs("species", local, attrs) == format_species("CFC-12")

    local = {"units": "ppt"} 
    attrs = {"units": "BLAH"}
    assert lookup_locals_and_attrs("units", local, attrs) == format_units("ppt")

    units = "ppt"
    attrs = {"units": "BLAH"}
    assert lookup_locals_and_attrs("units", locals(), attrs) == format_units("ppt")

    calibration_scale = "SIO-05"
    attrs = {"calibration_scale": "BLAH"}
    assert lookup_locals_and_attrs("calibration_scale", locals(), attrs) == format_calibration_scale("SIO-05")


def test_calibration_scale_default():

    # Read scale_defaults.csv file
    import pandas as pd

    with open_data_file("scale_defaults.csv", network="agage_test") as f:
        df = pd.read_csv(f, index_col="Species")
    
    # Remove whitespace
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)

    # Test species name
    assert calibration_scale_default("agage_test", "CH4") == df.loc["ch4", "calibration_scale"]
    assert calibration_scale_default("agage_test", "CH3Cl") == df.loc["ch3cl", "calibration_scale"]