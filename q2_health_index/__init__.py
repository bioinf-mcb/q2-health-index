from ._gmhi import calculate_gmhi
from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

__all__ = ['calculate_gmhi']