import pandas as pd
import numpy as np
import networkx as nx

from agage_archive import open_data_file
from agage_archive.data_selection import calibration_scale_default
from agage_archive.formatting import format_species


def scale_graph(species):
    """Build an directed graph from scale_convert.csv

    Args:
        species (str): Species
        
    Returns:
        nx.Graph: Undirected graph
    """

    with open_data_file("scale_convert.csv") as f:
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
        scale_new = calibration_scale_default(ds.attrs["network"], format_species(species))

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
