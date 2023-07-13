
from agage_archive.processing import lookup_locals_and_attrs, format_species, format_units, format_calibration_scale


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

