import plotly.graph_objects as go

from agage_archive.definitions import instrument_type_definition, unit_translator


colours = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
colour_counter = 0
colour_max = len(colours)

instrument_number, instrument_number_string = instrument_type_definition()

def plot_add_trace(fig, ds,
                name="", mode="lines"):

    global colour_counter

    # If data density is more than one point every hour, thin the dataset
    time_diff = ds.time.diff(dim="time")

    # Calculate the average time difference in seconds
    avg_time_diff = time_diff.mean().values / 1e9  
    if avg_time_diff < 600:
        print("Thinning dataset")
        ds_plot = ds.isel(time=slice(None, None, 60))
    else:
        ds_plot = ds.copy()

    inlets = set(ds.inlet_height.values)

    for inlet in inlets:

        # Get the indices of the inlet
        ind = ds_plot.inlet_height == inlet

        # Add trace
        fig.add_trace(
            go.Scatter(
                visible=True,
                mode=mode,
                marker=dict(color=colours[colour_counter % colour_max], size=5, symbol='cross'),
                line=dict(color=colours[colour_counter % colour_max], width=2),
                name=name + f", {inlet} m",
                x=ds_plot.time[ind].values,
                y=ds_plot.mf[ind].values,)
            )

        colour_counter += 1

    return fig


def plot_combined(ds, fig, mode="lines"):
    """ Plot multiple instruments on the same plot
    
    Args:
        ds (xarray.Dataset): Dataset containing multiple instruments (instrument_type variable needed)
        fig (plotly.graph_objects.Figure): Figure object
        mode (str) : how to plot the data. Defaults to "line", but can be set to "markers" if desired

    Returns:
        plotly.graph_objects.Figure: Figure object
    """

    global colour_counter

    # Get unique instrument types
    instrument_types = set(ds.instrument_type.values)

    # Remove unidentified instrument type
    instrument_types.discard(-1)

    # Create a trace for each instrument type
    for instrument_type in sorted(list(instrument_types)):

        instrument_type_name = list(instrument_number.keys())[list(instrument_number.values()).index(instrument_type)]

        # Get the indices of the instrument type
        ind = ds.instrument_type == instrument_type

        # Add trace
        fig = plot_add_trace(fig, ds.isel(time=ind),
                             name=f"{ds.attrs['site_code']}, {instrument_type_name}", mode=mode)

        colour_counter += 1
        
    return fig


def plot_single(ds, fig, mode="lines"):
    """ Plot a single instrument on the same plot

    Args:
        ds (xarray.Dataset): Dataset containing a single instrument
        fig (plotly.graph_objects.Figure): Figure object
        mode (str) : how to plot the data. Defaults to "line", but can be set to "markers" if desired
        
    Returns:
        plotly.graph_objects.Figure: Figure object
    """

    global colour_counter

    # Add trace
    fig = plot_add_trace(fig, ds,
                         name=f"{ds.attrs['site_code']}", mode=mode)

    colour_counter += 1
    
    return fig


def plot_datasets(datasets, mode="lines"):
    """ Plot datasets

    Args:
        datasets (list): List of xarray datasets containing AGAGE mole fraction data

    Returns:
        plotly.graph_objects.Figure: Figure object
    """

    unit = list(unit_translator.keys())[list(unit_translator.values()).index(datasets[0].mf.units)]

    # Create figure
    fig = go.Figure()

    # Set y-axis title to be species and units
    fig.update_yaxes(title_text=f"{datasets[0].attrs['species']} ({unit})")

    # Make figure less wide
    fig.update_layout(width=600,
                    height=500)

    # Make margins smaller
    fig.update_layout(margin=dict(l=20, r=20, t=20, b=20))

    # Move legend to top-right corner inside plot area
    fig.update_layout(legend=dict(x=0.02, y=1, yanchor="top"))

    # Change theme to simple
    fig.layout.template = "simple_white"

    for ds in datasets:

        # If instrument_type variable is present, split by instrument_type
        if "instrument_type" in ds.variables:

            fig = plot_combined(ds, fig, mode=mode)

        else:
            
            fig = plot_single(ds, fig, mode=mode)

    # Make sure legend is visible
    fig.update_layout(showlegend=True)

    # Convert to figure widget
    fig = go.FigureWidget(fig)

    fig.update_layout(xaxis=dict(tickformat='%Y-%m-%d'))

    return fig
