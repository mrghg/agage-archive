from pathlib import Path as _Path
import tarfile
from zipfile import ZipFile
import yaml
from fnmatch import fnmatch


__version__ = "0.0.1"


class Paths():

    def __init__(self,
                network = "",
                this_repo = False,
                errors = "raise"):
        """Class to store paths to data folders
        
        Args:
            network (str, optional): Network name. Defaults to "".
            this_repo (bool, optional): If True, look for the root and data folder within this repository (no config).
                If False, will look for the root and data folders and config file in the working directory.
                Defaults to False.
            errors (str, optional): If "raise", raise FileNotFoundError if file not found. If "ignore", return path

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
                "Config file not found. Try running util.setup first")

        # Read config file
        with open(self.config_file) as f:
            config = yaml.safe_load(f)
        
        self.user = config["user"]["name"]

        # If network is not set, exit
        if not network:
            return

        # Read all sub-paths associated with network
        for key, value in config["paths"][network].items():
            self.__setattr__(key, value)
            # Test that path exists
            full_path = self.data / network / value
            if not (full_path).exists() and errors == "raise":
                raise FileNotFoundError(f"Folder or zip archive {full_path} doesn't exist")
            # Test that path is either a folder or a zip archive
            if not (full_path.is_dir() or full_path.suffix == ".zip") and errors == "raise":
                raise FileNotFoundError(f"{full_path} is not a folder or zip archive")


def data_file_list(network = "",
                   sub_path = "",
                   pattern = "*",
                   ignore_hidden = True,
                   errors="raise"):
    """List files in data directory. Structure is data/network/sub_path
    sub_path can be a zip archive

    Args:
        network (str, optional): Network. Defaults to "".
        sub_path (str, optional): Sub-path. Defaults to "".
        pattern (str, optional): Pattern to match. Defaults to "*".

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

    pth = data_file_path("", network=network, sub_path=sub_path, errors=errors)

    if pth.suffix == ".zip":
        
        # If zip archive doesn't exist, return empty list
        if not pth.exists() and errors == "ignore":
            return network, return_sub_path(pth), []
        
        with ZipFile(pth, "r") as z:
            files = []
            for f in z.filelist:
                if fnmatch(f.filename, pattern) and not (ignore_hidden and f.filename.startswith(".")):
                    files.append(f.filename)
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
        return network, return_sub_path(pth), files

    

def data_file_path(filename,
                   network = "",
                   sub_path = "",
                   this_repo = False,
                   errors = "raise"):
    """Get path to data file. Structure is data/network/sub_path
    sub_path can be a zip archive, in which case the path to the zip archive is returned

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

    paths = Paths(network, this_repo=this_repo, errors=errors)

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
            raise FileNotFoundError(f"Can't find {filename} in {pth}")
    else:
        return pth / filename


def open_data_file(filename,
                   network = "",
                   sub_path = "",
                   verbose = False,
                   this_repo = False):
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

    pth = data_file_path("", network=network, sub_path=sub_path, this_repo=this_repo)
    
    if verbose:
        print(f"... opening {pth / filename}")

    if pth.suffix == ".zip":
        with ZipFile(pth, "r") as z:
            return z.open(filename)
    elif "tar.gz" in filename:
        return tarfile.open(pth / filename, "r:gz")
    else:
        return (pth / filename).open("rb")
