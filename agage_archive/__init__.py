from pathlib import Path as _Path
import configparser
import yaml

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


# class Paths():
#     def __init__(self, test=False):
#         """Class to store paths to data folders
#         """

#         # Get repository root
#         self.root = get_path().parent
#         self.data = self.root / "data"

#         if test == False:

#             # Check if config file exists
#             config_file = get_path("config.ini")
#             if not config_file.exists():
#                 raise FileNotFoundError(
#                     "Config file not found. Try running util.setup first")

#             # Read config file
#             config = configparser.ConfigParser()
#             config.read(get_path("config.ini"))

#             self.agage_gcmd = _Path(config["Paths"]["agage_md_path"])
#             #self.agage_gcmd = self.agage / "data-nc"
#             self.agage_gcms = _Path(config["Paths"]["agage_gcms_path"])
#             self.ale = _Path(config["Paths"]["ale_path"])
#             self.gage = _Path(config["Paths"]["gage_path"])

#             self.output = _Path(config["Paths"]["output_path"])

#         else:

#             # Use test data
#             self.agage_gcmd = self.root / "data/agage_test/data-nc"
#             self.agage_gcms = self.root / "data/agage_test/data-gcms-nc"
#             self.ale = self.root / "data/ale_gage_sio1993/ale"
#             self.gage = self.root / "data/ale_gage_sio1993/gage"

#             self.output = self.root / "data/agage_test/output"

#         # Check that data folders are there
#         for pth in [self.agage_gcmd,
#                     self.agage_gcms,
#                     self.ale,
#                     self.gage]:

#             if not pth.exists():
#                 raise FileNotFoundError(
#                     f"Folder {pth} doesn't exist")
        
#         # Check that NetCDF folders contain .nc files
#         for pth in [self.agage_gcmd, self.agage_gcms]:
#             if pth.suffix != ".zip":
#                 if not list(pth.glob("*.nc")):
#                     raise FileNotFoundError(
#                         f"""{pth} directory doesn't
#                         contain any NetCDF files"""
#                     )

#         # Check that ALE and GAGE folder contain tar.gz files
#         for pth in [self.ale, self.gage]:
#             if not list(pth.glob("*.gtar.gz")):
#                 raise FileNotFoundError(
#                     f"""{pth} directory doesn't
#                     contain any GCWerks .gtar.gz files"""
#                 )
