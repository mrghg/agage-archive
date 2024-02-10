import pandas as pd
import numpy as np
import networkx as nx
import json

from agage_archive.config import open_data_file
from agage_archive.data_selection import calibration_scale_default
from agage_archive.formatting import format_species


def resample(ds,
            resample_period = 3600,
            resample_threshold = 600,
            ):
    """Resample the dataset to a regular time interval
    
    Args:
        ds (xarray.Dataset): Dataset
        resample_period (int, optional): Period to resample to, in seconds. Defaults to 3600.
        resample_threshold (int, optional): Threshold for resampling, in seconds. Defaults to 600.
            If the median time difference is greater than this threshold, the dataset is not resampled.
            So, if the threshold is 600s and the resample_period is 3600s and 1-minute data is provided, 
            the dataset is resampled to houry. If 20-minute data is provided, the dataset is not resampled.

    Returns:
        xarray.Dataset: Resampled dataset
    """

    # Read variables.json
    with open_data_file("variables.json", this_repo=True) as f:
        variable_defaults = json.load(f)

    # Add some additional variables that can be processed, but don't need to be in variables.json
    variable_defaults["baseline"] = {"resample_method": ""}

    # check if median time difference is less than minimum_averaging_period
    if (ds.time.diff("time").median() < np.timedelta64(resample_threshold, "s")):

        df = ds.to_dataframe()
        df_resample = df.resample(f"{resample_period}S")
        df_resample_means = df_resample.mean()

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
            else:
                if var == "baseline":
                    # If any value in a resampled period is not 1, set baseline to 0
                    ds[var].values = np.where(df_resample[var].prod() != 1, 0, 1)
                elif var == "sampling_period":
                    # Set sampling period to 10 minutes
                    ds[var].values = np.ones_like(ds.time.values).astype(float) * resample_period
                elif var == "time":
                    pass
                else:
                    raise ValueError(f"Resample method not defined for {var}")

    else:
        return ds


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
