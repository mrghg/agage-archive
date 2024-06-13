from pathlib import Path as _Path
import tarfile
from zipfile import ZipFile
import yaml
from fnmatch import fnmatch, filter
import psutil
from shutil import copy


class Paths():

    def __init__(self,
                network = "",
                this_repo = False,
                errors = "raise",
                public = True):
        """Class to store paths to data folders
        
        Args:
            network (str, optional): Network name. Defaults to "".
            this_repo (bool, optional): If True, look for the root and data folder within this repository (no config).
                If False, will look for the root and data folders and config file in the working directory.
                Defaults to False.
            errors (str, optional): If "raise", raise FileNotFoundError if file not found. 
                If "ignore", return path.
                If "ignore_inputs", ignore errors in input paths. 
                Defaults to "raise".
            public (bool, optional): If True, output path is taken from output_path.
                If False, output path is taken from output_path_private. Defaults to True.

        Raises:
            FileNotFoundError: If config file doesn't exist
            FileNotFoundError: If folder or zip archive doesn't exist
            FileNotFoundError: If folder or zip archive is not a folder or zip archive   
        """

        # Get repository root
        # Do this by finding the location of the .git folder in the working directory
        # and then going up one level
        if not this_repo:
            working_directory = _Path.cwd()
            while True:
                if (working_directory / ".git").exists():
                    break
                else:
                    working_directory = working_directory.parent
                    if working_directory == _Path("/"):
                        raise FileNotFoundError("Can't find repository root")
        else:
            working_directory = _Path(__file__).parent.parent

        # Within working directory find package folder 
        # by looking for folder name with "_archive" in it, and __init__.py
        for pth in working_directory.glob("*_archive"):
            if (pth / "__init__.py").exists():
                self.root = pth
                break
        else:
            raise FileNotFoundError("Can't find package folder. Make sure your package has '_archive' in the folder name and __init__.py")

        self.data = self.root.parent / "data"

        # If this_repo is set, exit, to avoid config file confusion
        if this_repo:
            return

        # Check if config file exists
        self.config_file = self.root / "config.yaml"
        if not self.config_file.exists():
            raise FileNotFoundError(
                "Config file not found. Try running config.setup first")

        # Read config file
        with open(self.config_file) as f:
            config = yaml.safe_load(f)
        
        self.user = config["user"]["name"]

        # If network is not set, exit
        if not network:
            return

        # Check if network exists in config file
        if not network in config["paths"].keys():
            if errors == "ignore":
                return
            else:
                raise KeyError(f"Network {network} not found in config file")

        # Read all sub-paths associated with network and INPUTS
        for key, value in config["paths"][network].items():

            # If key is either of the output paths, go to next key
            if "output_path" in key:
                continue

            self.__setattr__(key, value)

            # Test that path exists
            if errors == "raise" or errors == "ignore_outputs":
                full_path = self.data / network / value
                if not (full_path).exists():
                    raise FileNotFoundError(f"Folder or zip archive {full_path} doesn't exist")
                if not (full_path.is_dir() or full_path.suffix == ".zip"):
                    raise FileNotFoundError(f"{full_path} is not a folder or zip archive")

        # Don't need to do the remaining checks if errors is set to ignore_outputs
        if "output_path" not in config["paths"][network]:
            if errors == "raise" or errors == "ignore_inputs":
                raise KeyError(f"Output path not found in config file")
            else:
                return
            
        # Set OUTPUT path
        if public:
            self.output_path = config["paths"][network]["output_path"]
        else:
            self.output_path = config["paths"][network]["output_path_private"]
        
        # Test that output path exists
        if errors == "raise" or errors == "ignore_inputs":
            full_path = self.data / network / self.output_path
            if not (full_path).exists():
                raise FileNotFoundError(f"Folder or zip archive {full_path} doesn't exist")
            if not (full_path.is_dir() or full_path.suffix == ".zip"):
                raise FileNotFoundError(f"{full_path} is not a folder or zip archive")


def setup():
    """ Setup the config.yml file for the agage_archive package.
    """

    #TODO: this_repo could cause confusion in the case that someone uses this script
    # to set up a config file in another repository
    # However, setting it for now, so that it doesn't look for config file before it's created
    paths = Paths(this_repo=True)

    header = '''# Use this file to store configuration settings
# All paths are relative to the network subfolder in the data directory
# If you need to put data files elsewhere, you'll need to use symlinks
---
'''

    config = {}

    # User name
    usr = input("Name (press enter for system ID):") or ""
    config["user"] = {"name": usr}

    config["paths"] = {
        "agage_test":
            {
                "md_path": "data-nc",
                "gcms_path": "data-gcms-nc",
                "gcms_flask_path": "data-gcms-flask-nc",
                "ale_path": "ale",
                "gage_path": "gage",
                "output_path": "output",
                "output_path_private": "output-private"
            },
        "agage":
            {
                "md_path": "data-nc.zip",
                "gcms_path": "data-gcms-nc.zip",
                "gcms_flask_path": "data-gcms-flask-nc.zip",
                "ale_path": "ale_gage_sio1993/ale",
                "gage_path": "ale_gage_sio1993/gage",
                "output_path": "agage-public-archive.zip",
                "output_path_private" : "agage-private-archive.zip"
        }
    }

    with open(paths.root / 'config.yaml', 'w') as configfile:
        # Write header lines
        configfile.write(header)
        # Dump config dictionary as yaml
        yaml.dump(config, configfile,
                  default_flow_style=False,
                  sort_keys=False)
        
    print(f"Config file written to {paths.root / 'config.yaml'}")
    print("Config file has been populated with default sub-paths relative to data/network. " + \
          "If you want to move the data elsewhere, manually modify the sub-paths in the config file. " + \
          "If the files need to be outside of the data directory, use symlinks.")


def data_file_list(network = "",
                   sub_path = "",
                   pattern = "*",
                   ignore_hidden = True,
                   errors="raise",
                   sub_directories = True):
    """List files in data directory. Structure is data/network/sub_path
    sub_path can be a zip archive

    Args:
        network (str, optional): Network. Defaults to "".
        sub_path (str, optional): Sub-path. Defaults to "".
        pattern (str, optional): Pattern to match. Defaults to "*".
        ignore_hidden (bool, optional): Ignore hidden files. Defaults to True.
        errors (str, optional): See options in Paths class. Defaults to "raise".
        sub_directories (bool, optional): If False, will remove sub-directories. Defaults to True.

    Returns:
        tuple: Tuple containing network, sub-path and list of files
    """

    def return_sub_path(full_path):
        pth = ""
        for p in full_path.parts[::-1]:
            if (p == "data") or (p == network):
                break
            else:
                pth = p + "/" + pth
        return pth

    def remove_sub_directories(files):
        nslash = [f.count("/") for f in files]
        max_slash = min(nslash)
        return [f for f in files if f.count("/") == max_slash]

    pth = data_file_path("", network=network, sub_path=sub_path, errors=errors)

    if pth.suffix == ".zip":
        
        # If zip archive doesn't exist, return empty list
        if not pth.exists() and "ignore" in errors:
            return network, return_sub_path(pth), []
        
        with ZipFile(pth, "r") as z:
            files = []
            for f in z.filelist:
                if fnmatch(f.filename, pattern) and not (ignore_hidden and f.filename.startswith(".")):
                    files.append(f.filename)
            if not sub_directories:
                files = remove_sub_directories(files)
            return network, return_sub_path(pth), files
    else:
        files = []
        # This is written this way to give the same output as the zip archive
        for f in pth.glob("**/*"):
            if fnmatch(str(f), "*" + pattern) and not (ignore_hidden and f.name.startswith(".")):
                # Append everything in file path after pth
                if f.is_dir():
                    files.append(str(f.relative_to(pth)) + "/")
                else:
                    files.append(str(f.relative_to(pth)))
        if not sub_directories:
            files = remove_sub_directories(files)
        return network, return_sub_path(pth), files


def data_file_path(filename,
                   network = "",
                   sub_path = "",
                   this_repo = False,
                   errors = "raise"):
    """Get path to data file. Structure is data/network/sub_path
    sub_path can be a zip archive, in which case the path to the zip archive is returned

    Note that by default, this function is only for input data.
    Output data is handled by output_path (which is a wrapper around this function)
    
    Args:
        filename (str): Filename
        network (str, optional): Network. Defaults to "".
        sub_path (str, optional): Sub-path. Defaults to ""
        this_repo (bool, optional): If True, look for the root and data folder within this repository (no config).
            If False, will look for the root and data folders and config file in the working directory.
        errors (str, optional): If "raise", raise FileNotFoundError if file not found. If "ignore", return path

    Raises:
        FileNotFoundError: Can't find file

    Returns:
        pathlib.Path: Path to file
    """

    paths = Paths(network,
                  this_repo=this_repo,
                  errors=errors)

    if network:
        pth = paths.data / network
    else:
        pth = paths.data

    if sub_path:
        pth = pth / sub_path
    
    if not pth.exists() and errors == "raise":
        raise FileNotFoundError(f"Can't find path {pth}")

    if pth.suffix == ".zip":
        # If filename is empty, user is just asking to return completed directory
        if filename == "":
            return pth
        # Otherwise, check if filename is in zip archive
        with ZipFile(pth, "r") as z:
            for f in z.filelist:
                if f.filename == filename:
                    return pth
            if errors == "raise":
                raise FileNotFoundError(f"Can't find {filename} in {pth}")
    else:
        return pth / filename


def open_data_file(filename,
                   network = "",
                   sub_path = "",
                   verbose = False,
                   this_repo = False,
                   errors = "raise"):
    """Open data file. Structure is data/network/sub_path
    sub_path can be a zip archive

    Args:
        filename (str): Filename
        network (str, optional): Network. Defaults to "".
        sub_path (str, optional): Sub-path. Defaults to "". Can be a zip archive or directory
        verbose (bool, optional): Print verbose output. Defaults to False.
        this_repo (bool, optional): If True, look for the root and data folder within this repository (no config).
            If False, will look for the root and data folders and config file in the working directory.

    Raises:
        FileNotFoundError: Can't find file

    Returns:
        file: File object
    """

    pth = data_file_path("", network=network,
                         sub_path=sub_path,
                         this_repo=this_repo,
                         errors=errors)
    
    if verbose:
        print(f"... opening {pth / filename}")

    if pth.suffix == ".zip":
        with ZipFile(pth, "r") as z:
            return z.open(filter(z.namelist(), filename)[0])
    elif "tar.gz" in filename:
        return tarfile.open(pth / filename, "r:gz")
    else:
        return (pth / filename).open("rb")


def output_path(network, species, site, instrument,
                extra = "", version="", public=True,
                errors="raise", network_out = ""):
    '''Determine output path and filename

    Args:
        network (str): Network
        species (str): Species
        site (str): Site
        instrument (str): Instrument
        extra (str, optional): Extra string to add to filename. Defaults to "".
        version (str, optional): Version number. Defaults to "".
        public (bool, optional): Whether the dataset is for public release. Default to True.
        errors (str, optional): How to handle errors if path doesn't exist. Defaults to "raise".
        network_out (str, optional): Network for filename. Defaults to "".

    Raises:
        FileNotFoundError: Can't find output path

    Returns:
        pathlib.Path: Path to output directory
        str: Filename
    '''

    # Get paths. Ignore errors since outputs may not exist at this stage
    paths = Paths(network, public=public, errors="ignore_outputs")

    version_str = f"_{version.replace(' ','')}" if version else ""
    
    # Can tweak data_file_path to get the output path
    output_path = data_file_path("", network = network,
                                sub_path = paths.output_path,
                                errors=errors)
    
    # Create filename
    if network_out:
        network_str = network_out.lower()
    else:
        network_str = network.lower()

    if instrument:
        instrument_str = f"_{instrument}"
    else:
        instrument_str = ""

    filename = f"{network_str}{instrument_str}_{site.lower()}_{species}{extra}{version_str}.nc"

    return output_path, filename


def copy_to_archive(src_file, network, public = True):
    """Copy file to archive. Structure is data/network/sub_path
    sub_path can be a zip archive

    Args:
        src_file (str): Source file
        network (str, optional): Network. Defaults to "".
        public (bool, optional): If True, copy to public archive. If False, copy to private archive.
            Defaults to True.
    Raises:
        FileNotFoundError: Can't find file

    Returns:
        file: File object
    """

    archive_path, _ = output_path(network, "_", "_", "_",
                                public=public)

    if archive_path.suffix == ".zip":
        with ZipFile(archive_path, "a") as z:
            z.write(src_file, arcname=src_file.name)
    else:
        # Copy file into pth directory
        copy(src_file, archive_path / src_file.name)


def is_jupyterlab_session():
    """Check whether we are in a Jupyter-Lab session.
    Taken from:
    https://stackoverflow.com/questions/57173235/how-to-detect-whether-in-jupyter-notebook-or-lab
    """

    # inspect parent process for any signs of being a jupyter lab server
    parent = psutil.Process().parent()
    if parent.name() == "jupyter-lab":
        return "jupyterlab"
    
    keys = (
        "JUPYTERHUB_API_KEY",
        "JPY_API_TOKEN",
        "JUPYTERHUB_API_TOKEN",
    )
    env = parent.environ()
    if any(k in env for k in keys):
        return "jupyterlab"

    return "notebook"


if __name__ == "__main__":
    setup()