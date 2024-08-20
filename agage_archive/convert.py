import pandas as pd
import numpy as np
import networkx as nx
import xarray as xr
import json
from openghg_calscales import convert as convert_scale

from agage_archive.config import open_data_file
from agage_archive.data_selection import calibration_scale_default
from agage_archive.formatting import format_species, format_variables, comment_append


def apply_resample_method(df, index, columns, variable_defaults, resample_period="3600s"):
    """Apply the resample method to the dataset

    Args:
        df (pd.DataFrame): DataFrame created from xarray dataset
        index (pd.DatetimeIndex): Index of the resampled DataFrame
        columns (list): List of columns in the DataFrame
        variable_defaults (dict): Variable defaults
        resample_period (str, optional): Period to resample to. Defaults to "3600s".

    Returns:
        pd.DataFrame: Resampled DataFrame
    """

    df_out = pd.DataFrame(index=index)

    for var in columns:
        if variable_defaults[var]["resample_method"] == "mean":
            df_out[var] = df[var].mean()
        elif variable_defaults[var]["resample_method"] == "median":
            df_out[var] = df[var].median()
        elif variable_defaults[var]["resample_method"] == "sum":
            df_out[var] = df[var].sum()
        elif variable_defaults[var]["resample_method"] == "standard_error":
            df_out[var] = df[var].mean() / np.sqrt(df[var].count())
        elif variable_defaults[var]["resample_method"] == "mode":
            if isinstance(df[var], pd.core.groupby.SeriesGroupBy):
                mode_value = df[var].apply(lambda x: x.mode().iloc[0] if not x.mode().empty else np.nan)
            else:
                mode_value = df[var].mode()
            df_out[var] = mode_value[0] if not mode_value.empty else np.nan
        else:
            if var == "baseline":
                # If any value in a resampled period is not 1, set baseline to 0
                df_out[var] = np.where(df[var].prod() != 1, 0, 1)
            elif var == "sampling_period":
                # Set sampling period to resample_period
                df_out[var] = pd.Timedelta(resample_period).total_seconds()
            elif var == "mf_variability":
                # Overwritten below
                pass
            else:
                raise ValueError(f"Resample method not defined for {var}")

    # Add in mole fraction standard deviation variable
    df_out["mf_variability"] = df["mf"].std().copy()

    if not "mf_count" in columns:
        df_out["mf_count"] = df["mf"].count().copy()
    # Else: Should have already been summed, as set in variables.json

    return df_out


def resampler(df, variable_defaults, resample_period="3600s"):
    """Resample the dataset to a regular time interval

    Args:
        df (pd.DataFrame): DataFrame created from xarray dataset
        variable_defaults (dict): Variable defaults
        resample_period (str, optional): Period to resample to. Defaults to "3600s". 
            Pandas alias for time period, e.g. "1H" for hourly, "1D" for daily.

    Returns:
        pd.DataFrame: Resampled DataFrame
    """

    df_resample = df.resample(resample_period, closed="left", label="left")

    df_out = apply_resample_method(df_resample,
                                   df_resample[df.columns[0]].mean().index,
                                   list(df.columns),
                                   variable_defaults,
                                   resample_period=resample_period)
    
    return df_out


def grouper(df, variable_defaults, resample_period="3600s"):
    """Group the dataset by inlet height or resample to a regular time interval, 
    if the inlet height changes more quickly than the resample period

    Args:
        df (pd.DataFrame): DataFrame created from xarray dataset
        variable_defaults (dict): Variable defaults
        resample_period (str, optional): Period to resample to. Defaults to "3600s". 
            Pandas alias for time period, e.g. "1H" for hourly, "1D" for daily.

    Returns:
        pd.DataFrame: Grouped/resampled DataFrame
    """

    # First, find the indices where the inlet height changes
    inlet_height_change_indices = np.where(np.concatenate([np.array([True]),
                                                           np.diff(df.inlet_height) != 0]))[0]
    inlet_height_change_indices = np.append(inlet_height_change_indices, len(df))

    dfs = []
    for i in range(len(inlet_height_change_indices)-1):

        df_slice = df.iloc[inlet_height_change_indices[i]:inlet_height_change_indices[i+1]]

        if np.unique(df_slice.inlet_height).size > 1:
            raise ValueError("Inlet height changes more than once in a single resample period. This must be a bug!")

        if (df_slice.index[-1] - df_slice.index[0]) > pd.Timedelta(resample_period):
            # Resample this slice
            df_avg = resampler(df_slice, variable_defaults, resample_period=resample_period)
            # If last resample period passes the end of the slice,
            # change the sampling period to the number of seconds between the last time point and the end of the slice
            if df_avg.index[-1] + pd.Timedelta(df_avg["sampling_period"].iloc[-1], unit="s") > \
                    df.index[inlet_height_change_indices[i+1]]:
                df_avg.loc[df_avg.index[-1], "sampling_period"] = (df.index[inlet_height_change_indices[i+1]] - \
                                                df_avg.index[-1]).seconds
        else:
            # Apply the resample method to the slice and store the result in a new dataframe
            df_avg = apply_resample_method(df_slice,
                                           pd.DatetimeIndex([df_slice.index[0]]),
                                           df_slice.columns,
                                           variable_defaults,
                                           resample_period=resample_period)
            if i == len(inlet_height_change_indices)-2:
                # If this is the last slice, set the sampling period to the number of seconds between the last time point
                # and the end of the dataset
                next_timestamp = df.index[-1] + pd.Timedelta(df["sampling_period"].iloc[-1], unit="s")
            else:
                next_timestamp = df.index[inlet_height_change_indices[i+1]]
            df_avg["sampling_period"] = (next_timestamp - df_slice.index[0]).seconds

        dfs.append(df_avg)

    df_avg = pd.concat(dfs, axis=0).sort_index()

    return df_avg


def resample(ds,
            resample_period = "3600s",
            resample_threshold = "600s",
            ignore_inlet = False
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
        
        # Remove duplicate time points, which seem to be there when inlets switch sometimes
        ds = ds.drop_duplicates("time", keep = "first")

        # Create new dataset to store resampled data
        ds_out = ds.isel(time=0)

        inlet_height_change_times = [ds.time.values[0]]
        inlet_height_change_times_delta = []

        for i in range(1, len(ds.time)):
            if ds.inlet_height.values[i] != ds.inlet_height.values[i-1]:
                inlet_height_change_times_delta.append(ds.time.values[i] - \
                                    inlet_height_change_times[-1])
                inlet_height_change_times.append(ds.time.values[i])

        inlet_height_change_times = pd.to_datetime(inlet_height_change_times)
        inlet_height_change_times_delta = pd.to_timedelta(inlet_height_change_times_delta)

        if inlet_height_change_times_delta.median() == pd.to_timedelta(0):

            df = resampler(ds.to_dataframe(), variable_defaults, resample_period=resample_period)
            comment_str = f"Resampled to {resample_period}."

        else:
            # Is inlet height changing more quickly than the resample period?
            if inlet_height_change_times_delta.median() < pd.to_timedelta(resample_period) and \
                    not ignore_inlet:

                df = grouper(ds.to_dataframe(), variable_defaults, resample_period=resample_period)
                comment_str = f"Grouped by inlet height and/or resampled to {resample_period}."

            else:

                df = resampler(ds.to_dataframe(), variable_defaults, resample_period=resample_period)
                comment_str = f"Resampled to {resample_period}."

        # Expand dimensions to include time
        # This step removes time attributes, so need to store and replace
        time_attrs = ds.time.attrs
        ds_out = ds_out.expand_dims(time=df.index)
        ds_out.time.attrs = time_attrs

        variables = list(ds.variables)
        variables.remove("time")

        for var in variables:
            ds_out[var] = xr.DataArray(df[var].values, dims=["time"],
                                    coords={"time": ds_out.time})
            ds_out[var].attrs = ds[var].attrs.copy()
        
        ds_out.time.attrs = ds.time.attrs.copy()

        ds_out["mf_variability"] = xr.DataArray(df["mf_variability"].values, dims=["time"],
                                            coords={"time": ds_out.time})
        ds_out["mf_variability"].attrs = variable_defaults["mf_variability"]["attrs"].copy()
        ds_out["mf_variability"].attrs["units"] = ds["mf"].attrs["units"]
        # Copy either "calibration_scale" or "scale" attribute from "mf" variable
        if "calibration_scale" in ds["mf"].attrs:
            ds_out["mf_variability"].attrs["calibration_scale"] = ds["mf"].attrs["calibration_scale"]
        elif "scale" in ds["mf"].attrs:
            ds_out["mf_variability"].attrs["calibration_scale"] = ds["mf"].attrs["scale"]
        else:
            raise ValueError("No calibration scale found for mole fraction variable.")

        ds_out["mf_count"] = xr.DataArray(df["mf_count"], dims=["time"],
                                coords={"time": ds_out.time})
        ds_out["mf_count"].attrs = variable_defaults["mf_count"]["attrs"].copy()
        ds_out["mf_count"].attrs["units"] = ""

        if "comment" not in ds.attrs:
            ds_out.attrs["comment"] = comment_str
        else:
            ds_out.attrs["comment"] = comment_append(ds_out.attrs["comment"], comment_str)
        if "comment" not in ds.time.attrs:
            ds_out.time.attrs["comment"] = comment_str
        else:
            ds_out.time.attrs["comment"] = comment_append(ds_out.time.attrs["comment"], comment_str)

        return ds_out

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
                              resample_threshold="1000D",
                              ignore_inlet=True) # Big number so that everything gets resampled

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
