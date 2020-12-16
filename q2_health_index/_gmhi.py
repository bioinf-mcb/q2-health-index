import pandas as pd


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

    # TODO CALCULATE GMHI (BELOW WE ARE CREATING FAKE DATA!!!)

    gmhi = pd.Series(data=[-3.438, -2.08342, 0.345],
                   index=['ACVD_1', 'ACVD_2', 'ACVD_3'], name='GMHI')

    gmhi_artifact = ctx.make_artifact('SampleData[AlphaDiversity]', gmhi)
    return gmhi_artifact
