from pathlib import Path as _pth

__version__ = "0.0.1"

_ROOT = _pth(__file__).parent


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
