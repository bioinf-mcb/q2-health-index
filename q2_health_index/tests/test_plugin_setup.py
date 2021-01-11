# -----------------------------------------------------------------------------
# Copyright (c) 2020-2021, Bioinformatics at Małopolska Centre of Biotechnology
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import unittest

from q2_health_index.plugin_setup import plugin as health_index_plugin


class PluginSetupTests(unittest.TestCase):

    def test_plugin_setup(self):
        self.assertEqual(health_index_plugin.name, 'health-index')
