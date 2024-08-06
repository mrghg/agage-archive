import pandas as pd
import numpy as np
import networkx as nx
import xarray as xr
import json
from openghg_calscales import convert as convert_scale

from agage_archive.config import open_data_file
from agage_archive.data_selection import calibration_scale_default
from agage_archive.formatting import format_species, format_variables


def resample(ds,
            resample_period = "3600s",
            resample_threshold = "600s",
            ):
    """Resample the dataset to a regular time interval
    
    Args:
        ds (xarray.Dataset): Dataset
        resample_period (str, optional): Period to resample to. Defaults to "3600s". 
            Pandas alias for time period, e.g. "1H" for hourly, "1D" for daily.
        resample_threshold (str, optional): Threshold for resampling, in seconds. Defaults to 600s.
            If the median time difference is greater than this threshold, the dataset is not resampled.
            So, if the threshold is 600s and the resample_period is 3600s and 1-minute data is provided, 
            the dataset is resampled to houry. If 20-minute data is provided, the dataset is not resampled.
            Pandas alias for time period, e.g. "600s" for 10 minutes.

    Returns:
        xarray.Dataset: Resampled dataset
    """

    # Read variables.json
    with open_data_file("variables.json", this_repo=True) as f:
        variable_defaults = json.load(f)

    # Read variable defaults for non-public data to find how to resample these variables
    with open_data_file("variables_not_public.json", this_repo=True) as f:
        variable_np_defaults = json.load(f)
    variable_defaults.update(variable_np_defaults)
    
    # check if median time difference is less than minimum_averaging_period
    if pd.to_timedelta(ds.time.diff("time").median().values) < \
            pd.to_timedelta(resample_threshold):

        # Pandas resampling is more efficient than xarray resampling, so use that
        df_resample = {}
        df_resample_means = {}
        df_resample_std = {}
        if "inlet" in ds.dims:
            inlets = ds.inlet.values
            for inlet in inlets:
                df = ds.sel(inlet=inlet).to_dataframe()
                df_resample[inlet] = df.resample(resample_period, closed="left", label="left")
                df_resample_means[inlet] = df_resample[inlet].mean()
                df_resample_std[inlet] = df_resample[inlet].std(ddof=0)
        else:
            inlets = [0]
            df = ds.to_dataframe()
            df_resample[0] = df.resample(resample_period, closed="left", label="left")
            df_resample_means[0] = df_resample[0].mean()
            df_resample_std[0] = df_resample[0].std(ddof=0)

        # Create list of variables from dataset, excluding time and inlet
        variables = list(ds.variables)
        variables.remove("time")
        if "inlet" in variables:
            variables.remove("inlet")

        # Create new dataset to store resampled data
        ds = ds.isel(dict(time=0)).copy()

        # Expand dimensions to include time
        # This step removes time attributes, so need to store and replace
        time_attrs = ds.time.attrs
        ds = ds.expand_dims(time=df_resample_means[inlets[0]].index).copy(deep=True)
        ds.time.attrs = time_attrs

        # Update variables with resampled data
        for var in variables:
            for inlet in inlets:
                if "inlet" in ds.dims:
                    inlet_dict = {"inlet": inlet}
                else:
                    inlet_dict = {}
                if variable_defaults[var]["resample_method"] == "mean":
                    ds[var].loc[inlet_dict] = df_resample_means[inlet][var].values.copy()
                elif variable_defaults[var]["resample_method"] == "median":
                    ds[var].loc[inlet_dict] = df_resample[inlet][var].median().values
                elif variable_defaults[var]["resample_method"] == "sum":
                    ds[var].loc[inlet_dict] = df_resample[inlet][var].sum().values
                elif variable_defaults[var]["resample_method"] == "standard_error":
                    ds[var].loc[inlet_dict] = df_resample_means[inlet][var].values / np.sqrt(df_resample[inlet][var].count().values)
                    ds[var].attrs["long_name"] = "Standard error in mean of " + ds[var].attrs["long_name"].lower()
                elif variable_defaults[var]["resample_method"] == "mode":
                    ds[var].loc[inlet_dict] = df_resample[inlet][var].apply(lambda x: x.mode()[0] if not x.mode().empty else -1).values
                else:
                    if var == "baseline":
                        # If any value in a resampled period is not 1, set baseline to 0
                        ds[var].loc[inlet_dict] = np.where(df_resample[inlet][var].prod() != 1, 0, 1)
                    elif var == "sampling_period":
                        # Set sampling period to resample_period
                        ds[var].loc[inlet_dict] = np.ones_like(ds.time.values).astype(float) * pd.Timedelta(resample_period).total_seconds()
                    elif var == "mf_variability":
                        # Overwritten below
                        pass
                    else:
                        raise ValueError(f"Resample method not defined for {var}")

        # Construct mf_variability dataArray to add to dataset
        mf_variability = np.array([df_resample_std[inlet]["mf"].values for inlet in inlets])

        # Add in mole fraction standard deviation variable
        ds["mf_variability"] = xr.DataArray(mf_variability.reshape(ds["mf"].shape),
                                            dims=ds.mf.dims,
                                            coords=ds.mf.coords)
        ds["mf_variability"].attrs = variable_defaults["mf_variability"]["attrs"].copy()
        ds["mf_variability"].attrs["units"] = ds["mf"].attrs["units"]
        # Copy either "calibration_scale" or "scale" attribute from "mf" variable
        if "calibration_scale" in ds["mf"].attrs:
            ds["mf_variability"].attrs["calibration_scale"] = ds["mf"].attrs["calibration_scale"]
        elif "scale" in ds["mf"].attrs:
            ds["mf_variability"].attrs["calibration_scale"] = ds["mf"].attrs["scale"]
        else:
            raise ValueError("No calibration scale found for mole fraction variable.")

        # Add in number of samples variable
        # If this variable is already in the dataset, it should have been resampled apprioriately
        if not "mf_count" in variables:
            mf_count = np.array([df_resample[inlet]["mf"].count().values for inlet in inlets])
            ds["mf_count"] = xr.DataArray(mf_count.reshape(ds["mf"].shape),
                                        dims=ds.mf.dims)
            ds["mf_count"].attrs = variable_defaults["mf_count"]["attrs"].copy()
            ds["mf_count"].attrs["units"] = ""
        
        # Append to comments saying that the data has been resampled
        if "comment" not in ds.attrs:
            ds.attrs["comment"] = f"Resampled to {resample_period}."
        else:
            ds.attrs["comment"] = ds.attrs["comment"] + f" Resampled to {resample_period}."
        if "comment" not in ds.time.attrs:
            ds.time.attrs["comment"] = f"Resampled to {resample_period}."
        else:
            ds.time.attrs["comment"] = ds.time.attrs["comment"] + f" Resampled to {resample_period}."

        return ds

    else:

        return ds


def monthly_baseline(ds, ds_baseline):
    '''Calculate monthly baseline mole fractions

    Args:
        ds (xr.Dataset): Dataset
        ds_baseline (xr.Dataset): Baseline dataset

    Returns:
        xr.Dataset: Dataset with monthly baseline mole fractions
    '''

    # Select baseline points
    ds_baseline_points = ds.where(ds_baseline.baseline == 1, drop=True)

    # Remove any baseline points where the mole fraction is NaN
    ds_baseline_points = ds_baseline_points.where(~np.isnan(ds_baseline_points.mf), drop=True)

    # Calculate monthly mean
    # If there are no baseline points, return an empty dataset with same attributes
    if len(ds_baseline_points.time) == 0:
        ds_monthly = ds.isel(time=[])
        ds_monthly["mf_variability"] = ds.mf.isel(time=[])
        ds_monthly["mf_variability"].attrs["units"] = ds.mf.attrs["units"]
        ds_monthly["mf_repeatability"] = ds.mf.isel(time=[])
        ds_monthly["mf_repeatability"].attrs["units"] = ds.mf.attrs["units"]
    else:
        ds_monthly = resample(ds_baseline_points,
                              resample_period="1MS",
                              resample_threshold="1000D") # Big number so that everything gets resampled

    # Relabel some attributes
    ds_monthly["mf_variability"].attrs["long_name"] = "Monthly standard deviation of baseline mole fractions"
    ds_monthly["mf_repeatability"].attrs["long_name"] = "Monthly standard error in mean of baseline mole fractions"
    
    # Copy global attributes
    ds_monthly.attrs = ds.attrs.copy()

    # Add baseline flag
    ds_monthly.attrs["baseline_flag"] = ds_baseline.attrs["baseline_flag"]

    # Run variable formatting
    ds_monthly = format_variables(ds_monthly,
                                  attribute_override={"mf_variability": ds_monthly["mf_variability"].attrs,
                                                      "mf_repeatability": ds_monthly["mf_repeatability"].attrs})

    return ds_monthly


def scale_convert(ds, scale_new):
    """Convert mole fraction from one scale to another

    Args:
        ds (xarray.Dataset): Dataset containing mole fractions
        scale_new (str): New scale to convert to. If None, no conversion is applied.
            If "default", the default scale for the species is used.
        
    Returns:
        ndarray,float: Mole fraction in new scale
    """

    # If no conversion required, return original dataset
    if scale_new == None:
        return ds

    # Find species
    species = ds.attrs["species"]

    # Get default scale, if needed
    if scale_new == "default":
        scale_new = calibration_scale_default(ds.attrs["network"],
                                              format_species(species))

    # If scales are the same, return original dataset
    if ds.attrs["calibration_scale"] == scale_new:
        return ds
    
    scale_original = ds.attrs["calibration_scale"]
    
    # Make a deep copy of the dataset
    ds_out = ds.copy(deep=True)

    # Apply conversion factor
    ds_out.mf.values = convert_scale(ds.mf, species, scale_original, scale_new).values

    # Update attributes
    ds_out.attrs["calibration_scale"] = scale_new
    for var in ds_out.variables:
        if "calibration_scale" in ds_out[var].attrs:
            ds_out[var].attrs["calibration_scale"] = scale_new
        if "scale" in ds_out[var].attrs:
            ds_out[var].attrs["scale"] = scale_new

    return ds_out


def separate_inlets(ds):
    """Separate dataset into dimension inlets based on inlet_height

    Args:
        ds (xarray.Dataset): Dataset

    Returns:
        xarray.Dataset: Dataset with inlet dimension
    """

    # If already separated by inlet, return dataset
    if "inlet" in ds.dims:
        return ds

    inlets = sorted(np.unique(ds.inlet_height.values))

    # Assign values to the inlet dimension based on inlet_height and replace all other values with NaN or equivalent
    ds_inlets = []
    for inlet in inlets:
        ds_inlet = ds.where(ds.inlet_height == inlet, drop=True)
        ds_inlet = ds_inlet.assign_coords({"inlet": inlet})
        ds_inlets.append(ds_inlet)

    ds = xr.concat(ds_inlets, dim="inlet")
    ds.inlet.attrs = {"long_name": "inlet height", "units": "m"}

    # If all values in one inlet are NaN, drop the inlet
    ds = ds.dropna(dim="inlet", how="all")

    # Remove inlet_height variable
    ds = ds.drop_vars("inlet_height")

    return ds
