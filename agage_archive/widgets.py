import xarray as xr
from IPython.display import clear_output

from agage_archive.config import Paths, data_file_list, open_data_file, is_jupyterlab_session
from agage_archive.visualise import plot_datasets


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
    options = sorted([f"{s}, {i}" for (s, i) in zip(sites, instruments) if "*" not in i])
    options += sorted([f"{s}, {i}" for (s, i) in zip(sites, instruments) if "*" in i])
    if instrument_site_dropdown:
        instrument_site_dropdown.options = filter_list(options)
    else:
        return filter_list(options)


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

    filenames = []
    for instrument_site in instrument_sites:
        site, instrument = instrument_site.split(', ')
        if "*" in instrument:
            filenames.append(f"{frequency}/{species}/individual/{instrument.split('*')[-1]}_{site}_{species}{frequency_suffix}_*.nc")
        else:
            filenames.append(f"{frequency}/{species}/{instrument}_{site}_{species}{frequency_suffix}_*.nc")

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

    filenames = get_filenames(species, frequency, instrument_site)
    
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

