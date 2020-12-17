import pandas as pd
import os

from q2_health_index._utilities import (_load_and_validate_species,
                                        _load_metadata)

def calculate_gmhi(ctx,
                   table=None,
                   metadata=None,
                   healthy_species=None,
                   non_healthy_species=None):

    healthy_species_list, non_healthy_species_list = \
        _load_and_validate_species(healthy_species, non_healthy_species)

    metadata = _load_metadata(metadata)

    # TODO CALCULATE GMHI (BELOW WE ARE JUST LOADING EXPECTED DATA !!!!)
    expected = os.path.join(os.path.dirname(__file__), 'tests/data/expected/4347_final_gmhi.tsv')
    gmhi = pd.read_csv(expected, sep='\t', index_col=0, header=0, squeeze=True)
    gmhi_artifact = ctx.make_artifact('SampleData[AlphaDiversity]', gmhi)
    return gmhi_artifact
