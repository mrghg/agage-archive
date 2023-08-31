from pathlib import Path as _Path
import tarfile
from zipfile import ZipFile
import yaml
from fnmatch import fnmatch


__version__ = "0.0.1"

_ROOT = _Path(__file__).parent


def get_path(sub_path=""):
    """Get path to data files

    Args:
        sub_path (str, optional): Sub-path. Defaults to "".

    Raises:
        Exception: If sub-path begins with '/'.

    Returns:
        pathlib.Path: Path to data files    
    """

    if sub_path:
        if sub_path[0] == "/":
            raise Exception("sub-path can't begin with '/'")

    path = _ROOT / sub_path

    return path


class Paths():

    def __init__(self, network = ""):
        """Class to store paths to data folders
        
        Args:
            network (str, optional): Network name. Defaults to "".

        Raises:
            FileNotFoundError: If config file doesn't exist
            FileNotFoundError: If folder or zip archive doesn't exist
            FileNotFoundError: If folder or zip archive is not a folder or zip archive   
        """

        # Get repository root
        self.root = get_path().parent
        self.data = self.root / "data"

        # Check if config file exists
        self.config_file = get_path("config.yaml")
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
            if not (full_path).exists():
                raise FileNotFoundError(f"Folder or zip archive {full_path} doesn't exist")
            # Test that path is either a folder or a zip archive
            if not (full_path.is_dir() or full_path.suffix == ".zip"):
                raise FileNotFoundError(f"{full_path} is not a folder or zip archive")


def data_file_list(network = "",
                   sub_path = "",
                   pattern = "*",
                   ignore_hidden = True):
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

    pth = data_file_path("", network=network, sub_path=sub_path)

    if pth.suffix == ".zip":
        with ZipFile(pth, "r") as z:
            files = []
            for f in z.filelist:
                if fnmatch(f.filename, pattern) and not (ignore_hidden and f.filename.startswith(".")):
                    files.append(f.filename)
            return network, return_sub_path(pth), files
    else:
        files = [f.name for f in pth.glob(pattern) if not (ignore_hidden and f.name.startswith("."))]
        return network, return_sub_path(pth), files
    

def data_file_path(filename,
                   network = "",
                   sub_path = ""):
    """Get path to data file. Structure is data/network/sub_path
    sub_path can be a zip archive, in which case the path to the zip archive is returned

    Args:
        filename (str): Filename
        network (str, optional): Network. Defaults to "".
        sub_path (str, optional): Sub-path. Defaults to "".

    Raises:
        FileNotFoundError: Can't find file

    Returns:
        pathlib.Path: Path to file
    """

    paths = Paths(network)

    if network:
        pth = paths.data / network
    else:
        pth = paths.data

    if sub_path:
        pth = pth / sub_path
    
    if not pth.exists():
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
                   verbose = False):
    """Open data file. Structure is data/network/sub_path
    sub_path can be a zip archive

    Args:
        filename (str): Filename
        network (str, optional): Network. Defaults to "".
        sub_path (str, optional): Sub-path. Defaults to "". Can be a zip archive or directory
        verbose (bool, optional): Print verbose output. Defaults to False.

    Raises:
        FileNotFoundError: Can't find file

    Returns:
        file: File object
    """

    pth = data_file_path("", network=network, sub_path=sub_path)
    
    if verbose:
        print(f"... opening {pth / filename}")

    if pth.suffix == ".zip":
        with ZipFile(pth, "r") as z:
            return z.open(filename)
    elif "tar.gz" in filename:
        return tarfile.open(pth / filename, "r:gz")
    else:
        return (pth / filename).open("rb")
