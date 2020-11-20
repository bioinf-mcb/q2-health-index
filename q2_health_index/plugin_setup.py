import q2_health_index

from q2_health_index._hello import print_hello
from qiime2.plugin import (Plugin, Citations)


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

plugin.visualizers.register_function(
    function=print_hello,
    inputs={},
    parameters={},
    input_descriptions={},
    parameter_descriptions={},
    name='Print hello.',
    description='Print hello on the screen.'
)
