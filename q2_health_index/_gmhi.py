import pandas as pd
import os
import biom

from q2_types.feature_table import (FeatureTable, Frequency)
from q2_health_index._utilities import (_load_and_validate_species,
                                        _load_metadata,
                                        _validate_metadata_is_superset,
                                        _validate_metadata_healthy_column)

def calculate_gmhi(ctx,
                   table=None,
                   metadata=None,
                   healthy_column=None,
                   healthy_states=None,
                   non_healthy_states=None,
                   healthy_species=None,
                   non_healthy_species=None):

    # load and validate species lists
    healthy_species_list, non_healthy_species_list = \
        _load_and_validate_species(healthy_species, non_healthy_species)

    # load and validate metadata
    metadata = _load_metadata(metadata)
    _validate_metadata_is_superset(metadata, table.view(biom.Table))

    # validate (non) healthy states and healthy column
    _validate_metadata_healthy_column(metadata, healthy_column,
                                      healthy_states, non_healthy_states)

    # TODO Pawel: move to separate function and add test (dada2_table.qza)
    # load and convert feature table (if needed)
    if table.type == FeatureTable[Frequency]:
        get_relative = ctx.get_action('feature_table', 'relative_frequency')
        table, = get_relative(table=table)

    # TODO Valentyn: CALCULATE GMHI (BELOW WE ARE JUST LOADING EXPECTED DATA !!!!)
    expected = os.path.join(os.path.dirname(__file__), 'tests/data/expected/4347_final_gmhi.tsv')
    gmhi = pd.read_csv(expected, sep='\t', index_col=0, header=0, squeeze=True)


    # TODO Pawel: add visualization

    gmhi_artifact = ctx.make_artifact('SampleData[AlphaDiversity]', gmhi)
    return gmhi_artifact
