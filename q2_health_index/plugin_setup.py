import q2_health_index

from q2_health_index._gmhi import calculate_gmhi
from qiime2.plugin import (Str, Plugin, Citations, Metadata, Visualization)
from q2_types.feature_table import (Frequency, FeatureTable)
from q2_types.feature_data import (Taxonomy, FeatureData)
from q2_types.sample_data import (SampleData, AlphaDiversity)

plugin = Plugin(
    name='health-index',
    version=q2_health_index.__version__,
    website="https://github.com/bioinf-mcb/q2-health-index",
    package='q2_health_index',
    citations=Citations.load('citations.bib', package='q2_health_index'),
    description=('This QIIME 2 plugin calculates the Gut Microbiome Health '
                 'Index (GMHI) created by Gupta et al. 2020 as well as '
                 'health and non-health prevalent taxa that can be later use '
                 'for calculating e.g. updated Gut Microbiome Health '
                 'Index .'),
    short_description='Plugin for calculating the Gut Microbiome Health Index (GMHI).'
)

plugin.pipelines.register_function(
    function=calculate_gmhi,
    inputs={'table': FeatureTable[Frequency]},
    parameters={
        'metadata': Metadata,
        'healthy_species': Str,
        'non_healthy_species': Str
    },
    outputs=[
        ('gmhi_results', SampleData[AlphaDiversity]),
        # ('gmhi_plot', Visualization)
    ],
    input_descriptions={'table': 'The feature frequency table to '
                                 'calculate Gut Microbiome Health Index from.'},
    parameter_descriptions={
        'metadata': 'Sample metadata containing individual_id_column '
                    'and healthy_column used in GMHI calculation.',
        'healthy_species': 'Path to file with healthy species.',
        'non_healthy_species': 'Path to file with non-healthy species.',
    },
    output_descriptions={
        'gmhi_results': 'Calculated GMHI in tabular form.',
        # 'gmhi_plot': 'Plot showing calculated GMHI distribution.'
    },
    name='Calculate GMHI',
    description='Calculate and plot Gut Microbial Health Index based on input data.'
)
