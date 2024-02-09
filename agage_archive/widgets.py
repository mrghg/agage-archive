import xarray as xr
import glob
from IPython.display import clear_output

from agage_archive.config import Paths, data_file_list, open_data_file
from agage_archive.visualise import plot_datasets


def file_search_species(network, species):
    """ Search for files containing species
    
    Args:
        species (str): Species to search for

    Returns:
        list: List of files containing species        
    """
    
    paths = Paths(network, errors="ignore_inputs")

    files = data_file_list(network,
                           sub_path=paths.output_path,
                           pattern=f"*/{species}/*.nc",
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


def update_instrument_site(change, network, instrument_site_dropdown):
    """ Update instrument and site dropdown

    Args:
        change (dict): Widget change dictionary
        network (str): Network
        network_site_dropdown (ipywidgets.Dropdown): Dropdown widget
    """

    files = file_search_species(network, change["new"])
    instruments, sites = instruments_sites(files)
    options = sorted([f"{s}, {i}" for (s, i) in zip(sites, instruments) if "*" not in i])
    options += sorted([f"{s}, {i}" for (s, i) in zip(sites, instruments) if "*" in i])
    if instrument_site_dropdown:
        instrument_site_dropdown.options = options
    else:
        return options


def get_filenames(species, instrument_sites):
    """ Get filenames from species and network/site
    
    Args:
        species (str): Species
        network_sites (list): List of network/site strings

    Returns:
        list: List of filenames
    """

    filenames = []
    for instrument_site in instrument_sites:
        site, instrument = instrument_site.split(', ')
        if "*" in instrument:
            filenames.append(f"*/{species}/individual/{instrument.split('*')[-1]}_{site}_{species}_*.nc")
        else:
            filenames.append(f"*/{species}/{instrument}_{site}_{species}_*.nc")

    return filenames


def load_datasets(network, filenames):
    """ Load datasets from filenames

    Args:
        filenames (list): List of filenames

    Returns:
        list: List of datasets
    """

    paths = Paths(network, errors="ignore_inputs")

    datasets = []
    for filename in filenames:
        with open_data_file(filename, network, paths.output_path, errors="ignore_inputs") as f:
            with xr.open_dataset(f) as ds:
                ds_species = ds.load()

        datasets.append(ds_species)

    return datasets


def plot_to_output(sender, network, species, network_site, output_widget):
    """ Plot to output widget

    Args:
        sender (ipywidgets.Button): Button widget
        species (str): Species
        network_site (str): Network and site
        output_widget (ipywidgets.Output): Output widget
    """

    if not network_site:
        with output_widget:
            clear_output(True)
            print("Please select a network and site") 

    filenames = get_filenames(species, network_site)
    datasets = load_datasets(network, filenames)

    with output_widget:
        clear_output()
        print(f"Plotting {species} for {network_site}... please wait...")
        clear_output(True)
        fig = plot_datasets(datasets)
        fig.show(renderer="notebook")


def show_netcdf_info(sender, network, species, network_site, output_widget):

    if not network_site:
        with output_widget:
            clear_output(True)
            print("Please select a network and site") 

    filenames = get_filenames(species, network_site)
    datasets = load_datasets(network, filenames)

    with output_widget:
        clear_output()
        for filename, dataset in zip(filenames, datasets):
            print(filename)
            print(dataset)
            print("-----------------------------------------")
            print("")

