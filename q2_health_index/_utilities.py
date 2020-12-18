import pandas as pd
import biom
import pkg_resources

from qiime2.plugin import Metadata

HEALTHY_SPECIES_DEFAULT_FP = \
    pkg_resources.resource_filename('q2_health_index', 'data/MH_species.txt')
NON_HEALTHY_SPECIES_DEFAULT_FP = \
    pkg_resources.resource_filename('q2_health_index', 'data/MN_species.txt')


def _load_and_validate_species(healthy_species_fp: str = None,
                               non_healthy_species_fp: str = None):
    healthy_species = _load_file(healthy_species_fp) \
        if healthy_species_fp else _load_file(HEALTHY_SPECIES_DEFAULT_FP)
    non_healthy_species = _load_file(non_healthy_species_fp) \
        if non_healthy_species_fp else _load_file(NON_HEALTHY_SPECIES_DEFAULT_FP)

    if not healthy_species:
        raise ValueError('Healthy species list is empty!')
    if not non_healthy_species:
        raise ValueError('Non-healthy species list is empty!')

    return healthy_species, non_healthy_species


def _load_file(file: str = None):
    with open(file, 'r') as f:
        return list(map(lambda x: x.strip(), f.readlines()))


def _load_metadata(metadata: Metadata = None):
    if not metadata:
        raise ValueError('Metadata parameter not provided!')
    metadata = metadata.to_dataframe()
    return metadata


# borrowed from q2_longitudinal
def _validate_metadata_is_superset(metadata: pd.DataFrame = None,
                                   table: biom.Table = None):
    metadata_ids = set(metadata.index.tolist())
    table_ids = set(table.ids())
    missing_ids = table_ids.difference(metadata_ids)
    if len(missing_ids) > 0:
        raise ValueError(f'Missing samples in metadata: {missing_ids}')


def _validate_and_extract_healthy_states(metadata: pd.DataFrame = None,
                                         healthy_column: str = None,
                                         healthy_states: str = None,
                                         non_healthy_states: str = None,
                                         ):
    # Basic validations
    if not healthy_column:
        raise ValueError('healthy_column parameter not provided!')
    if not healthy_states:
        raise ValueError('healthy_states parameter not provided!')
    if not non_healthy_states:
        raise ValueError('non_healthy_states parameter not provided!')
    if healthy_column not in metadata.columns:
        raise ValueError(f'\'{healthy_column}\' is not a column in your metadata.')
    # States validation
    if not healthy_states == 'rest':
        healthy_states = list(set(healthy_states.split(",")))
        for state in healthy_states:
            if state not in metadata[healthy_column].values:
                raise ValueError(f'Healthy state \'{state}\' is not represented '
                                 f'by any members of \'{healthy_column}\' column '
                                 'in metadata. Consider using a different '
                                 'healthy_column or state value.')
    if not non_healthy_states == 'rest':
        non_healthy_states = list(set(non_healthy_states.split(",")))
        for state in non_healthy_states:
            if state not in metadata[healthy_column].values:
                raise ValueError(f'Non-healthy state \'{state}\' is not represented '
                                 f'by any members of \'{healthy_column}\' column '
                                 'in metadata. Consider using a different '
                                 'healthy_column or state value.')
    if set(healthy_states) == set(non_healthy_states):
        raise ValueError(f'healthy_states and non_healthy_states '
                         f'parameters cannot be equal.')
    # TODO Consider a case when metadata dataset is superior wrt feature table
    # and contains additional (not relevant) state values
    if not healthy_states == 'rest' and not non_healthy_states == 'rest':
        number_of_state_values = sum([(metadata[healthy_column] == i).sum()
                                      for i in healthy_states + non_healthy_states])
        if number_of_state_values != len(metadata):
            raise ValueError(f'Number of healthy and non-healthy state values '
                             f'is not equal to the number of rows in metadata.')
    return healthy_states, non_healthy_states
