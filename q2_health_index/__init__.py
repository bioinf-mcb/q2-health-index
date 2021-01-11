# -----------------------------------------------------------------------------
# Copyright (c) 2020-2021, Bioinformatics at Małopolska Centre of Biotechnology
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from ._gmhi import calculate_gmhi
from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

__all__ = ['calculate_gmhi']
