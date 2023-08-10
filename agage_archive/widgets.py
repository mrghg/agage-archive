import xarray as xr
from IPython.display import clear_output, display, HTML
import plotly.offline as pyo

from agage_archive import Paths
from agage_archive.visualise import plot_datasets

paths = Paths()


def file_search_species(species):
    """ Search for files containing species
    
    Args:
        species (str): Species to search for

    Returns:
        list: List of files containing species        
    """
    species_path = paths.output / species
    files = species_path.glob('*.nc')
    return list(files)


def networks_sites(files):
    """ Get networks and sites from files

    Args:
        files (list): List of files

    Returns:
        tuple: Tuple of lists containing networks and sites
    """

    networks = []
    sites = []
    for file in files:
        networks.append(file.stem.split('_')[0])
        sites.append(file.stem.split('_')[1])
    return networks, sites


def update_network_site(change, network_site_dropdown):
    """ Update network and site dropdown

    Args:
        change (dict): Change dictionary
        network_site_dropdown (ipywidgets.Dropdown): Dropdown widget
    """

    files = file_search_species(change["new"])
    networks, sites = networks_sites(files)
    network_site_dropdown.options = sorted([f"{s}, {n}" for (s, n) in zip(sites, networks)])


def get_filenames(species, network_sites):
    """ Get filenames from species and network/site
    
    Args:
        species (str): Species
        network_sites (list): List of network/site strings

    Returns:
        list: List of filenames
    """

    filenames = []
    for network_site in network_sites:
        site, network = network_site.split(', ')
        filenames.append(f"{species}/{network}_{site}_{species}.nc")

    return filenames


def load_datasets(filenames):
    """ Load datasets from filenames

    Args:
        filenames (list): List of filenames

    Returns:
        list: List of datasets
    """

    datasets = []
    for filename in filenames:
        datasets.append(xr.open_dataset(paths.output / filename))

    return datasets


def plot_to_output(sender, species, network_site, output_widget):
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
            #display(output_widget)

    filenames = get_filenames(species, network_site)
    datasets = load_datasets(filenames)

    with output_widget:
        clear_output()
        print(f"Plotting {species} for {network_site}... please wait...")
        clear_output(True)
        fig = plot_datasets(datasets)
        fig.show(renderer="notebook")
