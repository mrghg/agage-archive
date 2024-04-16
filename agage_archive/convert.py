import pandas as pd
import numpy as np
import networkx as nx
import xarray as xr
import json

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

        # Pandas does resampling more efficiently, for some reason
        df = ds.to_dataframe()
        df_resample = df.resample(resample_period, closed="left", label="left")
        df_resample_means = df_resample.mean()
        df_resample_std = df_resample.std(ddof=0) # Use biased estimator for standard deviation

        # Create new dataset to store resampled data
        ds = ds.isel(time=0)
        ds = ds.expand_dims(time=df_resample_means.index)

        # Create list of variables from dataset, excluding time
        variables = list(ds.variables)
        variables.remove("time")

        # Update variables with resampled data
        for var in variables:
            if variable_defaults[var]["resample_method"] == "mean":
                ds[var].values = df_resample_means[var].values
            elif variable_defaults[var]["resample_method"] == "median":
                ds[var].values = df_resample[var].median().values
            elif variable_defaults[var]["resample_method"] == "sum":
                ds[var].values = df_resample[var].sum().values
            elif variable_defaults[var]["resample_method"] == "standard_error":
                ds[var].values = df_resample[var].mean().values / np.sqrt(df_resample[var].count().values)
                ds[var].attrs["long_name"] = "Standard error in mean of " + ds[var].attrs["long_name"].lower()
            elif variable_defaults[var]["resample_method"] == "mode":
                ds[var].values = df_resample[var].apply(lambda x: x.mode()[0] if not x.mode().empty else -1)
            else:
                if var == "baseline":
                    # If any value in a resampled period is not 1, set baseline to 0
                    ds[var].values = np.where(df_resample[var].prod() != 1, 0, 1)
                elif var == "sampling_period":
                    # Set sampling period to resample_period
                    ds[var].values = np.ones_like(ds.time.values).astype(float) * pd.Timedelta(resample_period).total_seconds()
                elif var == "mf_variability":
                    # Overwritten below
                    pass
                else:
                    raise ValueError(f"Resample method not defined for {var}")

        # Add in mole fraction standard deviation variable
        ds["mf_variability"] = xr.DataArray(df_resample_std["mf"].values, dims=["time"],
                                            coords={"time": ds.time})
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
        if not "mf_N" in variables:
            ds["mf_N"] = xr.DataArray(df_resample["mf"].count().values, dims=["time"],
                                    coords={"time": ds.time})
            ds["mf_N"].attrs = variable_defaults["mf_N"]["attrs"].copy()
            ds["mf_N"].attrs["units"] = ""
        
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

    return ds_monthly


def scale_graph(species):
    """Build an directed graph from scale_convert.csv

    Args:
        species (str): Species
        
    Returns:
        nx.Graph: Undirected graph
    """

    with open_data_file("scale_convert.csv", this_repo=True) as f:
        data = pd.read_csv(f)

    # Check if species is in the data
    if species not in data['Species'].values:
        raise ValueError(f"Species {species} not found in scale_convert.csv")

    # Filter the data for the specified species
    species_data = data[data['Species'] == species].iloc[0]

    G = nx.DiGraph()

    # Add edges to the graph based on the species data
    for col in data.columns[1:]:
        source, target = col.split('/')
        if not pd.isna(species_data[col]):
            G.add_edge(target, source, weight=species_data[col])
            G.add_edge(source, target, weight=1/species_data[col])

    return G


def scale_convert(ds, scale_new):
    """Convert mole fraction from one scale to another

    Args:
        ds (xarray.Dataset): Dataset containing mole fractions
        scale_new (str): New scale to convert to. If None, no conversion is applied.
            If "default", the default scale for the species is used.
        
    Returns:
        ndarray,float: Mole fraction in new scale
    """

    def n2o_scale_function(time, invert = False):
        """Function to apply to N2O mole fractions to convert from SIO-93 to SIO-98

        Args:
            time (pd.Timestamp): Timestamp
            mf (ndarray): Mole fractions
            invert (bool, optional): If True, apply the inverse function. Defaults to False.

        Returns:
            ndarray: Mole fractions adjusted for time-variation between scales (excluding factor)
        """

        # calculate days elapsed since 2nd March 1978
        days_since_ale_start = (time - pd.Timestamp("1978-03-02")).days.values

        a0=1.00664
        a1=-0.56994e-3
        a2=-0.65398e-3
        a3=0.13083e-3
        a4=-0.20742e-4

        t = (days_since_ale_start-3227.)/365.25
        f = 1./(a0 + a1*t + a2*t**2 + a3*t**3 + a4*t**4)
        if invert:
            f = 1./f

        # Apply f to mf only between 1st May 1984 and 31st March 1990
        f_out = np.ones_like(f).astype(float)
        idx = (days_since_ale_start>=2252) & (days_since_ale_start<=4412)
        f_out[idx] = f[idx]

        return f_out

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
    
    # TODO: CHECK IF SCALES ARE IN SCALE_CONVERT.CSV

    # Make a deep copy of the dataset
    ds_out = ds.copy(deep=True)

    G = scale_graph(format_species(species))

    # Find the path and calculate the conversion factor
    try:
        path = nx.shortest_path(G, source=scale_original,
                                target=scale_new,
                                weight='weight')
        conversion_factor = np.ones_like(ds.time.values).astype(float)
        for i in range(len(path)-1):
            if format_species(species) == "n2o":
                if path[i] == "SIO-93" and path[i+1] == "SIO-98":
                    conversion_factor *= n2o_scale_function(ds.time.to_series().index)
                elif path[i] == "SIO-98" and path[i+1] == "SIO-93":
                    conversion_factor *= n2o_scale_function(ds.time.to_series().index, invert=True)
            conversion_factor *= G[path[i]][path[i+1]]['weight']

    except nx.NetworkXNoPath:
        raise ValueError(f"No conversion path found between {scale_original} and {scale_new} for {species}.")

    # Apply conversion factor
    ds_out.mf.values *= conversion_factor
    ds_out.mf_repeatability.values *= conversion_factor

    # Update attributes
    ds_out.attrs["calibration_scale"] = scale_new

    return ds_out
