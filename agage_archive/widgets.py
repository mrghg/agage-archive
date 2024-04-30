import xarray as xr
from IPython.display import clear_output, display
import ipywidgets as widgets
from glob import glob
from pathlib import Path

from agage_archive.config import Paths, data_file_list, open_data_file, is_jupyterlab_session
from agage_archive.visualise import plot_datasets


# Global dictionary to store filenames associated with instrument/site string displayed in dropdown
instrument_site_filenames = {}


def file_search_species(network, frequency, species, public = True):
    """ Search for files containing species
    
    Args:
        network (str): Network
        frequency (str): Frequency ("event" or "monthly")
        species (str): Species to search for
        public (bool): Search public or private archive

    Returns:
        list: List of files containing species        
    """
    
    paths = Paths(network, errors="ignore_inputs")

    if public:
        output_path = paths.output_path
    else:
        output_path = paths.output_path_private

    files = data_file_list(network,
                           sub_path=f"{output_path}",
                           pattern=f"{frequency}/{species}/*.nc",
                           errors="ignore_inputs")[2]

    # Remove baseline files
    files = [f for f in files if "baseline" not in f]

    return sorted(files)


def instruments_sites(files):
    """ Get networks and sites from files

    Args:
        files (list): List of files

    Returns:
        tuple: Tuple of lists containing networks and sites
    """

    instruments = []
    sites = []
    individual = []

    for file in files:
        if "individual" in file:
            individual = "*"
        else:
            individual = ""
        filename = file.split('/')[-1]
        instruments.append(individual + filename.split('_')[0])
        sites.append(filename.split('_')[1])

    return instruments, sites


def update_instrument_site(species,
                        frequency,
                        network,
                        public,
                        instrument_site_dropdown):
    """ Update instrument and site dropdown

    Args:
        species (str): Species
        frequency (str): Frequency ("event" or "monthly")
        network (str): Network
        public (str): Load from public or private archive (public or private)
        instrument_site_dropdown (ipywidgets.Dropdown): Dropdown widget
    """
    def find_duplicates(options):
        """ Find duplicates in options list
        """
        seen = set()
        duplicates = set()
        for option in options:
            if option in seen:
                duplicates.add(option)
            else:
                seen.add(option)
        
        error_str = "Duplicate options found in instrument and site combinations:"

        if duplicates:
            for duplicate in duplicates:
                error_str += f"\n{duplicate}"
            raise ValueError(error_str)

    global instrument_site_filenames

    # Clear contents of instrument_site_filenames
    instrument_site_filenames.clear()

    files = file_search_species(network, frequency, species,
                                public = {"public": True, "private": False}[public])
    instruments, sites = instruments_sites(files)

    options = []
    options_files = []
    options_individual = []
    options_individual_files = []

    for f, i, s in zip(files, instruments, sites):
        if "*" in i:
            options_individual.append(f"{s}, {i}")
            options_individual_files.append(f)
        else:
            options.append(f"{s}, {i}")
            options_files.append(f)

    # Check for duplicates and raise an error if found
    find_duplicates(options)
    find_duplicates(options_individual)
    
    # Populate dictionary with instrument/site strings as keys and filenames as values
    for option in sorted(options):
        instrument_site_filenames[option] = options_files[options.index(option)]
    for option in sorted(options_individual):
        instrument_site_filenames[option] = options_individual_files[options_individual.index(option)]

    # Update dropdown widget if passed as kwarg, otherwise return options list
    if instrument_site_dropdown:
        instrument_site_dropdown.options = instrument_site_filenames.keys()
    else:
        return instrument_site_filenames.keys()


def get_filenames(frequency, instrument_sites):
    """ Get filenames from species and network/site
    
    Args:
        frequency (str): Frequency
        instrument_sites (list): List of instrument/site strings

    Returns:
        list: List of filenames
    """

    global instrument_site_filenames

    # Extract everything in file path from frequency onwards
    filenames = []
    for instrument_site in instrument_sites:
        filepath = instrument_site_filenames[instrument_site]
        filenames.append(f"{frequency}/{filepath.split(frequency + '/')[1]}")

    return filenames


def load_datasets(network, filenames, public = True):
    """ Load datasets from filenames

    Args:
        network (str): Network
        filenames (list): List of filenames
        public (bool): Load from public or private archive

    Returns:
        list: List of datasets
    """

    paths = Paths(network, errors="ignore_inputs")

    if public:
        output_path = paths.output_path
    else:
        output_path = paths.output_path_private

    datasets = []
    for filename in filenames:
        with open_data_file(filename, network, output_path, errors="ignore_inputs") as f:
            with xr.open_dataset(f) as ds:
                ds_species = ds.load()

        datasets.append(ds_species)

    return datasets


def plot_to_output(sender, network, frequency, species, instrument_site, public,
                   output_widget):
    """ Plot to output widget

    Args:
        sender (ipywidgets.Button): Button widget
        network (str): Network
        frequency (str): Frequency
        species (str): Species
        network_site (str): Network and site
        output_widget (ipywidgets.Output): Output widget
    """

    if not instrument_site:
        with output_widget:
            clear_output(True)
            print("Please select a network and site") 

    filenames = get_filenames(frequency, instrument_site)

    datasets = load_datasets(network, filenames,
                            public = {"public": True, "private": False}[public])

    renderer = is_jupyterlab_session()

    with output_widget:
        clear_output()
        print(f"Plotting {species} for {instrument_site}... please wait...")
        clear_output(True)
        fig = plot_datasets(datasets)
        fig.show(renderer=renderer)


def show_netcdf_info(sender, network, frequency, instrument_site, public,
                    output_widget):
    """ Show NetCDF info to output widget

    Args:
        sender (ipywidgets.Button): Button widget
        species (str): Species
        frequency (str): Frequency ("event" or "monthly")
        network_site (str): Network and site
        public (str): Load from public or private archive (public or private)
        output_widget (ipywidgets.Output): Output widget
    """

    if not instrument_site:
        with output_widget:
            clear_output(True)
            print("Please select a network and site") 

    filenames = get_filenames(frequency, instrument_site)
    datasets = load_datasets(network, filenames,
                            public = {"public": True, "private": False}[public])

    with output_widget:
        clear_output()
        for filename, dataset in zip(filenames, datasets):
            print(filename)
            print(dataset)
            print("-----------------------------------------")
            print("")


def dashboard(network,
            frequencies = ["event", "monthly"]):
    """ Create dashboard for visualising data

    Args:
        network (str): Network
        frequencies (list): List of frequencies ("event" or "monthly")
    """
    
    paths = Paths(network, errors="ignore_inputs")    
    
    # Get species names from the output directory structure
    species = []
    for f in data_file_list(network, paths.output_path, errors="ignore_inputs")[2]:
        if "/" in f:
            species.append(f.split("/")[1])
        else:
            continue

    if len(species) == 0:
        raise ValueError("No files found in the output directory")

    species = sorted(set(species))

    # Public or private archive radio button
    public_button = widgets.RadioButtons(
        options=["public", "private"],
        description='Archive:',
        disabled=False,
        default="public"
    )

    # Create dropdown widget
    species_dropdown = widgets.Dropdown(
        options=species,
        description='Species:',
        disabled=False,
        default=species[0]
        )

    # Create file_type widget
    frequency_dropdown = widgets.Dropdown(
        options=frequencies,
        description='Frequency:',
        disabled=False,
        default=frequencies[0]
    )

    # Selection widget for network and site
    instrument_site = widgets.SelectMultiple(
        options=update_instrument_site(species[0],
                                    frequencies[0],
                                    network,
                                    "public",
                                    None),
        description='Site, instrument:',
        disabled=False,
        indent=True,
        style={'description_width': 'initial'}
    )

    # Plotting button
    plot_button = widgets.Button(description="Plot")

    # Output widget
    output = widgets.Output()
    output_netcdf = widgets.Output()

    # Text widget to explain what the asterisk means
    asterisk_text = widgets.HTML(value="<p>* Asterisk indicates individual, rather than combined file</p>")

    # Update network and site dropdown when species is changed
    species_dropdown.observe(lambda change:
                            update_instrument_site(change["new"],
                                                frequency_dropdown.value,
                                                network,
                                                public_button.value,
                                                instrument_site),
                            names="value")

    frequency_dropdown.observe(lambda change:
                            update_instrument_site(species_dropdown.value,
                                            change["new"],
                                            network,
                                            public_button.value,
                                            instrument_site),
                            names="value")

    public_button.observe(lambda change:
                        update_instrument_site(species_dropdown.value,
                                        frequency_dropdown.value,
                                        network,
                                        change["new"],
                                        instrument_site),
                        names="value")

    # Plot to output when button is clicked
    plot_button.on_click(lambda x: plot_to_output(x, network,
                                                frequency_dropdown.value,
                                                species_dropdown.value,
                                                instrument_site.value,
                                                public_button.value,
                                                output))
    plot_button.on_click(lambda x: show_netcdf_info(x, network,
                                                    frequency_dropdown.value,
                                                    instrument_site.value,
                                                    public_button.value,
                                                    output_netcdf))

    display(public_button)
    display(species_dropdown)
    display(frequency_dropdown)
    display(instrument_site)
    display(asterisk_text)
    display(plot_button)
    display(output)
    display(output_netcdf)