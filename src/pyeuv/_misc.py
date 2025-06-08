import functools
import xarray as xr
from importlib_resources import files


@functools.cache
def _read_coeffs(file):
    return xr.open_dataset(files('pyeuv._coeffs').joinpath(file))


def get_euv91_coeffs():
    return (_read_coeffs('euv91_bands_dataset.nc').copy(), _read_coeffs('euv91_lines_dataset.nc').copy(),
            _read_coeffs('euv91_full_dataset.nc').copy())
