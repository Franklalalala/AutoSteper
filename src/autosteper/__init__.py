import sys

if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = "AutoSteper"
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

from .optimizers import ASE_Optimizer, Gaussian_Optimizer, XTB_Optimizer, Multi_Optimizer
from .Autosteper import AutoSteper

__all__ = ['ASE_Optimizer', 'Gaussian_Optimizer', 'XTB_Optimizer', 'Multi_Optimizer', 'AutoSteper']
