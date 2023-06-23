from agage_archive import io
import pandas as pd

paths = io.Paths()

def test_scale_convert():

    scale_conversion = pd.read_csv(paths.root / "data/scale_convert.csv",
                                   index_col="Species")

    mf = 1.

    # Test a set of conversion factors
    mf_new = io.scale_convert("CFC-11", "SIO-93", "SIO-98", mf)
    assert mf_new / mf == scale_conversion.loc["CFC-11", "SIO-98/SIO-93"]

    mf_new = io.scale_convert("CFC-11", "SIO-93", "SIO-05", mf)
    assert mf_new / mf == scale_conversion.loc["CFC-11", "SIO-98/SIO-93"] * \
        scale_conversion.loc["CFC-11", "SIO-05/SIO-98"]

    mf_new = io.scale_convert("CCl4", "SIO-93", "SIO-05", mf)
    assert mf_new / mf == scale_conversion.loc["CCl4", "SIO-98/SIO-93"] * \
        scale_conversion.loc["CCl4", "SIO-05/SIO-98"]

    mf_new = io.scale_convert("CH4", "CSIRO-94", "TU-87", mf)
    assert mf_new / mf == scale_conversion.loc["CH4", "TU-87/CSIRO-94"]
    
