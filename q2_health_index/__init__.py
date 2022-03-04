# -----------------------------------------------------------------------------
# Copyright (c) 2020-2021, Bioinformatics at Małopolska Centre of Biotechnology
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from ._gmhi import gmhi_fit, gmhi_predict, gmhi_predict_viz
from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

__all__ = ['gmhi_fit', 'gmhi_predict', 'gmhi_predict_viz']
