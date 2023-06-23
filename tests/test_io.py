from agage_archive import io
import pandas as pd

paths = io.Paths()

def test_scale_convert():

    scale_conversion = pd.read_csv(paths.root / "data/scale_convert.csv",
                                   index_col="Species")

    t = pd.Timestamp("1991-01-01")
    mf = 1.

    # Test a set of conversion factors
    mf_new = io.scale_convert("CFC-11", "SIO-93", "SIO-98", t, mf)
    assert mf_new / mf == scale_conversion.loc["CFC-11", "SIO-98/SIO-93"]

    mf_new = io.scale_convert("CFC-11", "SIO-93", "SIO-05", t, mf)
    assert mf_new / mf == scale_conversion.loc["CFC-11", "SIO-98/SIO-93"] * \
        scale_conversion.loc["CFC-11", "SIO-05/SIO-98"]

    mf_new = io.scale_convert("CCl4", "SIO-93", "SIO-05", t, mf)
    assert mf_new / mf == scale_conversion.loc["CCl4", "SIO-98/SIO-93"] * \
        scale_conversion.loc["CCl4", "SIO-05/SIO-98"]

    mf_new = io.scale_convert("CH4", "CSIRO-94", "TU-87", t, mf)
    assert mf_new / mf == scale_conversion.loc["CH4", "TU-87/CSIRO-94"]

    mf_new = io.scale_convert("N2O", "SIO-93", "SIO-05", t, mf)
    assert mf_new[0] / mf == scale_conversion.loc["N2O", "SIO-98/SIO-93"] * \
        scale_conversion.loc["N2O", "SIO-05/SIO-98"]
    
    # Period where N2O time conversion applies
    mf_new = io.scale_convert("N2O", "SIO-93", "SIO-05", pd.Timestamp("1985-01-01"), mf)
    assert mf_new[0] / mf == scale_conversion.loc["N2O", "SIO-98/SIO-93"] * \
        0.9962230167482587

