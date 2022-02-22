# -----------------------------------------------------------------------------
# Copyright (c) 2020-2021, Bioinformatics at Małopolska Centre of Biotechnology
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import q2_health_index

from q2_health_index._gmhi import gmhi_predict, gmhi_predict_viz
from qiime2.plugin import (Int, Str, Float, Plugin, Citations, Metadata,
                           Visualization)
from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_types.sample_data import SampleData, AlphaDiversity


basic_parameters = {
        'healthy_species_fp': Str,
        'non_healthy_species_fp': Str,
        'mh_prime': Int,
        'mn_prime': Int,
        'rel_thresh': Float,
        'log_thresh': Float,
    }

basic_parameters_descriptions = {
        'healthy_species_fp': 'Path to file with healthy species (taxonomy ' 
                              'is based on MetaPhlAn 2).',
        'non_healthy_species_fp': 'Path to file with non-healthy species ' 
                                  '(taxonomy is based on MetaPhlAn 2).',
        'mh_prime': 'Median from the top 1% healthy samples in training  '
                    'dataset (see Gupta et al. 2020 Methods section).',
        'mn_prime': 'Median from the top 1% non-healthy samples in training '
                    'dataset (see Gupta et al. 2020 Methods section).',
        'rel_thresh': 'Relative frequency based threshold for discarding '
                      'insignificant OTU.',
        'log_thresh': 'Normalization value for log10 in the last step of '
                      'GMHI calculation.',
    }

plugin = Plugin(
    name='health-index',
    version=q2_health_index.__version__,
    website="https://github.com/bioinf-mcb/q2-health-index",
    package='q2_health_index',
    citations=Citations.load('citations.bib', package='q2_health_index'),
    description=('This QIIME 2 plugin calculates and visualizes the Gut '
                 'Microbiome Health Index (GMHI) according to the algorithm '
                 'from Gupta et al. 2020.'),
    short_description='Plugin for calculating the Gut Microbiome Health '
                      'Index (GMHI).'
)

plugin.pipelines.register_function(
    function=gmhi_predict,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency]},
    parameters=basic_parameters,
    outputs=[
        ('gmhi_results', SampleData[AlphaDiversity]),
    ],
    input_descriptions={'table': 'The feature frequency table to calculate '
                                 'Gut Microbiome Health Index from.'},
    parameter_descriptions=basic_parameters_descriptions,
    output_descriptions={
        'gmhi_results': 'Calculated GMHI in tabular form.',
    },
    name='Calculate GMHI',
    description='Calculate Gut Microbial Health Index based on input data. '
)

plugin.pipelines.register_function(
    function=gmhi_predict_viz,
    inputs={
        'table': FeatureTable[Frequency | RelativeFrequency],
    },
    parameters={
        **basic_parameters,
        'metadata': Metadata,
    },
    outputs=[
        ('gmhi_results', SampleData[AlphaDiversity]),
        ('gmhi_plot', Visualization)
    ],
    input_descriptions={'table': 'The feature frequency table to calculate '
                                 'Gut Microbiome Health Index from.',
                        },
    parameter_descriptions={
        **basic_parameters_descriptions,
        'metadata': 'Metadata used for visualization [REQUIRED].',
    },
    output_descriptions={
        'gmhi_results': 'Calculated GMHI in tabular form.',
        'gmhi_plot': 'Bar plot showing calculated GMHI distribution.'
    },
    name='Calculate & visualize GMHI',
    description='Calculate and plot Gut Microbial Health Index based on '
                'input data and metadata. '
)
