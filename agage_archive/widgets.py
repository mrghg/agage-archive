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

    print(files)

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
    def filter_list(input_list):
        """Filter list to keep only first duplicate"""
        seen_set = set()
        result_list = []

        for item in input_list:
            if item not in seen_set:
                seen_set.add(item)
                result_list.append(item)

        return result_list

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

    if len(set(options)) != len(options):
        # Find and print duplicates
        seen = set()
        duplicates = set()
        for option in options:
            if option in seen:
                duplicates.add(option)
            else:
                seen.add(option)
        print(duplicates)
        
        raise ValueError("Duplicate options found in instrument and site combinations")
    if len(set(options_individual)) != len(options_individual):
        print(options_individual)
        raise ValueError("Duplicate options found in individual instrument and site combinations")

    for option in sorted(options):
        instrument_site_filenames[option] = options_files[options.index(option)]
    for option in sorted(options_individual):
        instrument_site_filenames[option] = options_individual_files[options_individual.index(option)]

    if instrument_site_dropdown:
        instrument_site_dropdown.options = instrument_site_filenames.keys()
    else:
        return instrument_site_filenames.keys()

    # options = sorted([f"{s}, {i}" for (s, i) in zip(sites, instruments) if "*" not in i])
    # options += sorted([f"{s}, {i}" for (s, i) in zip(sites, instruments) if "*" in i])
    # if instrument_site_dropdown:
    #     instrument_site_dropdown.options = filter_list(options)
    # else:
    #     return filter_list(options)


def get_filenames(species, frequency, instrument_sites):
    """ Get filenames from species and network/site
    
    Args:
        species (str): Species
        frequency (str): Frequency
        instrument_sites (list): List of instrument/site strings

    Returns:
        list: List of filenames
    """

    if frequency == "monthly":
        frequency_suffix = "-monthly"
    else:
        frequency_suffix = ""

    # filenames = []
    # for instrument_site in instrument_sites:
    #     site, instrument = instrument_site.split(', ')
    #     if "*" in instrument:
    #         search_str = f"{frequency}/{species}/individual/{instrument.split('*')[-1]}_{site}_{species}{frequency_suffix}_*.nc"
    #     else:
    #         search_str = f"{frequency}/{species}/{instrument}_{site}_{species}{frequency_suffix}_*.nc"

    #     files = glob(search_str)
    #     if len(files) == 0:
    #         raise ValueError(f"No files found for {instrument_site}")
    #     elif len(files) > 1:
    #         raise ValueError(f"Multiple files found for {instrument_site}")
    #     else:
    #         # Append just the filename without the path
    #         filenames.append(files[0].split('/')[-1])
    
    print(instrument_site_filenames)

    #if instrument_site_filenames:
    filenames = [Path(instrument_site_filenames[instrument_site]).name for instrument_site in instrument_sites]
    # else:
    #     filenames = []

    return filenames

# # Test for get_filenames
# def test_get_filenames():

#     species = "ch4"
#     frequency = "event"
#     instrument_sites = ["CGO, combined"]

#     filenames = get_filenames(species, frequency, instrument_sites)

#     assert len(filenames) == 1
#     assert filenames[0] == "AGAGE_mace_head_ch4-monthly_2000-01.nc"

#     print("get_filenames passed")

# test_get_filenames()


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

    filenames = get_filenames(species, frequency, instrument_site)
    print(filenames)
    datasets = load_datasets(network, filenames,
                            public = {"public": True, "private": False}[public])

    renderer = is_jupyterlab_session()

    with output_widget:
        clear_output()
        print(f"Plotting {species} for {instrument_site}... please wait...")
        clear_output(True)
        fig = plot_datasets(datasets)
        fig.show(renderer=renderer)


def show_netcdf_info(sender, network, frequency, species, instrument_site, public,
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

    filenames = get_filenames(species, frequency, instrument_site)
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
                                                    species_dropdown.value,
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