import configparser
from pathlib import Path
import xarray as xr

from agage_archive import get_path


class Paths():
    def __init__(self):

        # Check if config file exists
        config_file = get_path("config.ini")
        if not config_file.exists():
            raise FileNotFoundError(
                "Config file not found. Try running util.setup first")

        # Read config file
        config = configparser.ConfigParser()
        config.read(get_path("config.ini"))

        self.input_path = Path(config["Paths"]["input_path"])
        self.output_path = Path(config["Paths"]["output_path"])

        self.data_suffix = "data-nc"
        self.data_gcms_suffix = "data-gcms-nc"

        if not (self.input_path / self.data_suffix).exists():
            raise FileNotFoundError(
                f"Data folder must contain {self.data_suffix} folder")
        else:
            if not list((self.input_path / self.data_suffix).glob("*.nc")):
                raise FileNotFoundError(
                    f"{self.data_suffix} directory doesn't contain any NetCDF files"
                )

        if not (self.input_path / self.data_gcms_suffix).exists():
            raise FileNotFoundError(
                f"Data folder must contain {self.data_gcms_suffix} folder")
        else:
            if not list((self.input_path / self.data_gcms_suffix).glob("*.nc")):
                raise FileNotFoundError(
                    f"""{self.data_gcms_suffix} directory doesn't
                    contain any NetCDF files"""
                )


paths = Paths()


def read_nc(species, site, instrument):

    nc_file = paths.input_path / paths.data_suffix
    nc_file = nc_file / f"AGAGE-{instrument}_{site}_{species}.nc"

    if not nc_file.exists():
        raise FileNotFoundError(f"Can't find file {nc_file}")

    with xr.open_dataset(nc_file) as f:
        ds = f.load()

    return ds

