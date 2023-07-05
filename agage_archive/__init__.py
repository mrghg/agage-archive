from pathlib import Path
import configparser

__version__ = "0.0.1"

_ROOT = Path(__file__).parent


def get_path(sub_path=""):
    """Get path to data files

    Parameters
    ----------
    sub_path : str
        path to data files, relative to py12box/data directory

    Returns
    -------
    pathlib.Path
        pathlib Path to data folder/file

    """

    if sub_path:
        if sub_path[0] == "/":
            raise Exception("sub-path can't begin with '/'")

    path = _ROOT / sub_path

    return path


class Paths():
    def __init__(self):
        """Class to store paths to data folders
        """

        # Get repository root
        self.root = get_path().parent

        # Check if config file exists
        config_file = get_path("config.ini")
        if not config_file.exists():
            raise FileNotFoundError(
                "Config file not found. Try running util.setup first")

        # Read config file
        config = configparser.ConfigParser()
        config.read(get_path("config.ini"))

        self.agage = Path(config["Paths"]["agage_path"])
        self.agage_gcmd = self.agage / "data-nc"
        self.agage_gcms = self.agage / "data-gcms-nc"
        self.ale = Path(config["Paths"]["ale_path"])
        self.gage = Path(config["Paths"]["gage_path"])

        self.output = Path(config["Paths"]["output_path"])

        # Check that data folders are there
        for pth in [self.agage_gcmd,
                    self.agage_gcms,
                    self.ale,
                    self.gage]:

            if not pth.exists():
                raise FileNotFoundError(
                    f"Folder {pth} doesn't exist")
        
        # Check that NetCDF folders contain .nc files
        for pth in [self.agage_gcmd, self.agage_gcms]:
            if not list(pth.glob("*.nc")):
                raise FileNotFoundError(
                    f"""{pth} directory doesn't
                    contain any NetCDF files"""
                )

        # Check that ALE and GAGE folder contain .C files
        for pth in [self.ale, self.gage]:
            if not list(pth.glob("*.gtar.gz")):
                raise FileNotFoundError(
                    f"""{pth} directory doesn't
                    contain any GCWerks .gtar.gz files"""
                )

