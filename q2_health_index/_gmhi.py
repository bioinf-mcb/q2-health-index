# -----------------------------------------------------------------------------
# Copyright (c) 2020-2021, Bioinformatics at Małopolska Centre of Biotechnology
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import qiime2
import numpy as np
import pandas as pd

from q2_types.feature_table import FeatureTable, Frequency, RelativeFrequency
from q2_health_index._utilities import (_load_and_validate_species,
                                        _load_metadata,
                                        _validate_metadata_is_superset)


def gmhi_fit(train_df: pd.DataFrame, 
            train_labels: np.ndarray, 
            theta_f: float, 
            theta_d: float):
    '''
    Fits the function to the available training set with parameters theta_f and theta_d
    Parameters
    ----------
    train_df : pd.DataFrame
        Dataframe with index as samples and columns as features (species).
    train_labels : np.ndarray
        An array for binary categorical target variable (healthy and npn-healthy in the article) 
        encoded as 0 and 1.
    theta_f : float
        Fold difference between proportion of samples with relevant species present in it for a bianry state. 
        Defines relevant features for binary states. 
    theta_f : float
        Difference between proportion of samples with relevant species present in it for a bianry state.
        Defines relevant species for binary states. 
    '''
    # Extracting healthy and non-healthy cohorts form dataframe
    Healthy = train_df.iloc[train_labels]
    Nonhealthy = train_df.iloc[~train_labels]

    # Transposing similarly to article 
    Healthy = Healthy.T
    Nonhealthy = Nonhealthy.T
    
    # Calculating the proportion of samples with species present
    PH = (Healthy > 0).sum(axis=1) * 100 / Healthy.shape[1]
    PNH = (Nonhealthy > 0).sum(axis=1) * 100 / Nonhealthy.shape[1]
    
    # Deriving fold and normal differences
    PH_diff = (PH-PNH)
    PH_fold = (PH/PNH)
    PNH_fold = (PNH/PH)
    
    # Masking only species defined as important for health state by parameters
    H_signature = train_df.loc[:, (PH_fold >= theta_f) & (PH_diff >= theta_d)].T
    NH_signature = train_df.loc[:, (PNH_fold >= theta_f) &(PH_diff <= -theta_d)].T
    
    # Extracting lists of species
    healthy_species_list = list(H_signature.index)
    nonhealthy_species_list = list(NH_signature.index)
    
    # Deriving counts
    H_sig_count = (H_signature > 0).sum(axis=0)
    H_sig_count.name = "H_sig_count"
    NH_sig_count = (NH_signature > 0).sum(axis=0)
    NH_sig_count.name = "NH_sig_count"

    constant = pd.concat([H_sig_count, NH_sig_count], axis=1)
  
    # MH_prime
    HC1 = constant.sort_values(by = ["H_sig_count", "NH_sig_count"], ascending=[False, True])
    H_constant = HC1["H_sig_count"][:(int(Healthy.shape[1] / 100))].median()
    
    # MN_prime
    NHC1 = constant.sort_values(by = ["H_sig_count", "NH_sig_count"], ascending=[True, False])
    NH_constant = NHC1["NH_sig_count"][:(int(Nonhealthy.shape[1] / 100))].median() 
    
    return (healthy_species_list, 
            nonhealthy_species_list,
            H_constant, 
            NH_constant)


def calculate_gmhi(ctx,
                   table=None,
                   healthy_species_fp=None,
                   non_healthy_species_fp=None,
                   mh_prime=7,
                   mn_prime=31,
                   rel_thresh=0.00001,
                   log_thresh=0.00001):

    # Load and validate species lists
    healthy_species_list, non_healthy_species_list = \
        _load_and_validate_species(healthy_species_fp, non_healthy_species_fp)

    # Load and convert feature table (if needed)
    if table.type == FeatureTable[Frequency]:
        get_relative = ctx.get_action('feature_table', 'relative_frequency')
        table, = get_relative(table=table)
    assert table.type == FeatureTable[RelativeFrequency], \
        'Feature table not of the type \'RelativeFrequency\''

    # Keep columns as samples, rows as taxonomies
    table_df = table.view(pd.DataFrame).T

    # Consider only species from the full taxonomy
    table_df.index = table_df.index.str.split(';').str[-1].str.strip()

    # Remove unclassified and virus species - suitable both for 16S and
    # Metagenome Sequencing if valid taxonomy is provided
    na_species = table_df.index.str.contains('unclassified|virus', regex=True)
    species_profile_2 = table_df[~na_species]

    # Re-normalization of species' relative abundances after removing
    # unclassified and virus species
    species_profile_3 = species_profile_2.apply(lambda x: x / x.sum(), axis=0)
    species_profile_3[species_profile_3 < rel_thresh] = 0

    # Extracting Health-prevalent species
    MH_species = species_profile_3[
        species_profile_3.index.isin(healthy_species_list)]
    # Extracting Health-scarce species
    MN_species = species_profile_3[
        species_profile_3.index.isin(non_healthy_species_list)]

    assert not MH_species.empty, \
        "Could not find healthy species in the feature table."
    assert not MN_species.empty, \
        "Could not find non-healthy species in the feature table."

    # Shannon index + Alpha diversity
    MH_not_zero = MH_species[MH_species > 0]
    MN_not_zero = MN_species[MN_species > 0]
    MH_shannon = MH_not_zero.apply(lambda x: -np.sum(np.log(x) * x), axis=0)
    MN_shannon = MN_not_zero.apply(lambda x: -np.sum(np.log(x) * x), axis=0)

    # Richness of Health-prevalent species
    R_MH = MH_not_zero.count()
    # Richness of Health-scarce species
    R_MN = MN_not_zero.count()

    psi_MH = (R_MH / mh_prime) * MH_shannon
    psi_MN = (R_MN / mn_prime) * MN_shannon

    gmhi_df = np.log10((psi_MH + log_thresh) / (psi_MN + log_thresh))
    gmhi_df.name = 'GMHI'

    # Create and return artifact
    gmhi_artifact = ctx.make_artifact('SampleData[AlphaDiversity]', gmhi_df)

    return gmhi_artifact


def calculate_gmhi_viz(ctx,
                       table=None,
                       metadata=None,
                       healthy_species_fp=None,
                       non_healthy_species_fp=None,
                       mh_prime=7,
                       mn_prime=31,
                       rel_thresh=0.00001,
                       log_thresh=0.00001):

    # Calculate GMHI
    gmhi_artifact = calculate_gmhi(ctx, table, healthy_species_fp,
                                   non_healthy_species_fp, mh_prime, mn_prime,
                                   rel_thresh, log_thresh)

    # Load metadata
    metadata_df = _load_metadata(metadata)
    # Limit metadata to samples preset in the feature table
    table_df = table.view(pd.DataFrame).T
    metadata_df = _validate_metadata_is_superset(metadata_df, table_df.T)
    metadata = qiime2.Metadata(metadata_df)

    # Create visualization (box plots) similar to that from alpha-diversity
    get_alpha_diversity_plot = ctx.get_action('diversity',
                                              'alpha_group_significance')
    gmhi_viz = get_alpha_diversity_plot(alpha_diversity=gmhi_artifact,
                                        metadata=metadata)

    return gmhi_artifact, gmhi_viz[0]
    