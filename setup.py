# ----------------------------------------------------------------------------
# Copyright (c) 2020, Bioinformatics at Małopolska Centre of Biotechnology.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
import versioneer

setup(
    name="q2-health-index",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    author="Tomasz Kościółek",
    author_email="tomasz.kosciolek@uj.edu.pl",
    description="Calculate the Health Index from microbiome data",
    license='BSD-3-Clause',
    url="https://qiime2.org",
    entry_points={
        'qiime2.plugins':
        ['q2-health-index=q2_health_index.plugin_setup:plugin']
    },
    package_data={
        'q2_health_index.tests': ['data/*'],
        # TODO: remove 'tests/data/expected/*' after implementing calculate_gmhi !!!!
        'q2_health_index': ['data/*', 'tests/data/expected/*',
                            'assets/index.html', 'citations.bib']},
    zip_safe=False,
)