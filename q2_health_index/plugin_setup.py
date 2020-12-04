`import q2_health_index

from q2_health_index._hello import print_hello
from qiime2.plugin import (Plugin, Citations, Metadata, Visualization)
from q2_types.feature_table import (RelativeFrequency)
from q2_types.feature_data import (Taxonomy)


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

plugin.methods.register_function(
    function=calculate_gmhi,
    inputs={'table': FeatureTable[Frequency],
            'list_healthy': FeatureData[Taxonomy],
            'list_non_healthy': FeatureData[Taxonomy]},
    parameters={'metadata': Metadata},
    outputs=[('gmhi_plot', Visualization),
             ('gmhi_results', Metadata)],
    input_descriptions={'table': 'The feature frequency table to calculate Gut Microbiome Health Index from.',
                        'list_healthy' : 'List of healthy species imported as FeatureData[Taxonomy] artifact',
                        'list_non_healthy' : 'List of non-healthy species imported as FeatureData[Taxonomy] artifact'},
    parameter_descriptions={'metadata': 'Sample metadata containing  '
                                        'individual_id_column, and other metadata for use in '
                                        'Gut Microbiome Health Index Calculation.'},
    name='Calculate GMHI',
    description='Calculate and plot Gut Microbial Health Index based on input data.'
)
