# -----------------------------------------------------------------------------
# Copyright (c) 2020-2021, Bioinformatics at Małopolska Centre of Biotechnology
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import unittest
from warnings import filterwarnings

import pandas as pd
import pandas.testing as pdt
import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import health_index

from q2_health_index._utilities import (_load_file,
                                        _load_metadata,
                                        _load_and_validate_species,
                                        _validate_metadata_is_superset,
                                        HEALTHY_SPECIES_DEFAULT_FP,
                                        NON_HEALTHY_SPECIES_DEFAULT_FP)

filterwarnings("ignore", category=UserWarning)
filterwarnings("ignore", category=RuntimeWarning)

HEALTHY_SPECIES_DEFAULT = _load_file(HEALTHY_SPECIES_DEFAULT_FP)
NON_HEALTHY_SPECIES_DEFAULT = _load_file(NON_HEALTHY_SPECIES_DEFAULT_FP)


class TestUtilities(TestPluginBase):
    package = 'q2_health_index.tests'

    def test_species_load_default(self):
        healthy, non_healthy = _load_and_validate_species(None, None)
        self.assertIsNotNone(HEALTHY_SPECIES_DEFAULT)
        self.assertIsNotNone(NON_HEALTHY_SPECIES_DEFAULT)
        self.assertListEqual(healthy, HEALTHY_SPECIES_DEFAULT)
        self.assertListEqual(non_healthy, NON_HEALTHY_SPECIES_DEFAULT)

    def test_species_load_healthy_wrong(self):
        with self.assertRaisesRegex(FileNotFoundError, "No such file or "
                                                       "directory"):
            infile = self.get_data_path("input/species/do-not-exists.txt")
            _load_and_validate_species(healthy_species_fp=infile)

    def test_species_load_non_healthy_wrong(self):
        with self.assertRaisesRegex(FileNotFoundError, "No such file or "
                                                       "directory"):
            infile = self.get_data_path("input/species/do-not-exists.txt")
            _load_and_validate_species(non_healthy_species_fp=infile)

    def test_species_load_healthy_empty(self):
        with self.assertRaisesRegex(ValueError, "Healthy species list is "
                                                "empty!"):
            infile = self.get_data_path("input/species/empty_species.txt")
            _load_and_validate_species(healthy_species_fp=infile)

    def test_species_load_non_healthy_empty(self):
        with self.assertRaisesRegex(ValueError, "Non-healthy species list is "
                                                "empty!"):
            infile = self.get_data_path("input/species/empty_species.txt")
            _load_and_validate_species(non_healthy_species_fp=infile)

    def test_species_load_healthy_fake(self):
        infile = self.get_data_path("input/species/fake_MH_species.txt")
        healthy, non_healthy = \
            _load_and_validate_species(healthy_species_fp=infile)
        self.assertListEqual(healthy, ['s__fake_1', 's__fake_2'])
        self.assertListEqual(non_healthy, NON_HEALTHY_SPECIES_DEFAULT)

    def test_species_load_non_healthy_fake(self):
        infile = self.get_data_path("input/species/fake_MN_species.txt")
        healthy, non_healthy = \
            _load_and_validate_species(non_healthy_species_fp=infile)
        self.assertListEqual(non_healthy, ['s__fake_non_1', 's__fake_non_2',
                                           's__fake_non_3'])
        self.assertListEqual(healthy, HEALTHY_SPECIES_DEFAULT)

    def test_species_load_fake(self):
        infile1 = self.get_data_path("input/species/fake_MH_species.txt")
        infile2 = self.get_data_path("input/species/fake_MN_species.txt")
        healthy, non_healthy = _load_and_validate_species(infile1, infile2)
        self.assertListEqual(healthy, ['s__fake_1', 's__fake_2'])
        self.assertListEqual(non_healthy, ['s__fake_non_1', 's__fake_non_2',
                                           's__fake_non_3'])

    def test_metadata_not_provided(self):
        with self.assertRaisesRegex(ValueError, "Metadata parameter not "
                                                "provided!"):
            _load_metadata(None)

    def test_metadata_load_simple(self):
        infile = self.get_data_path("input/metadata/simple_metadata.tsv")
        metadata = _load_metadata(qiime2.Metadata.load(infile))
        metadata_exp = pd.read_csv(infile, sep='\t')
        self.assertListEqual(list(metadata.columns), ['Age', 'Healthy'])
        self.assertListEqual(list(metadata.index),
                             list(metadata_exp['sample-id']))
        self.assertListEqual(list(metadata.Healthy),
                             list(metadata_exp['Healthy']))
        self.assertListEqual(list(metadata.Age),
                             list(metadata_exp['Age']))

    def test_metadata_validate_simple(self):
        table_file = self.get_data_path("input/abundances"
                                        "/simple_relative_abundances.qza")
        metadata_file = self.get_data_path(
            'input/metadata/simple_metadata.tsv')
        metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
        table = qiime2.Artifact.load(table_file).view(pd.DataFrame)
        metadata_new = _validate_metadata_is_superset(metadata, table)
        self.assertListEqual(sorted(metadata_new.index), sorted(table.index))

    def test_metadata_validate_simple_wrong(self):
        with self.assertRaisesRegex(ValueError,
                                    "Missing samples in metadata: {'MOCK-"):
            table_file = self.get_data_path("input/abundances/dada2_table.qza")
            metadata_file = self.get_data_path("input/metadata"
                                               "/simple_metadata.tsv")
            metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
            table = qiime2.Artifact.load(table_file).view(pd.DataFrame)
            _validate_metadata_is_superset(metadata, table)


class TestPredictGmhi(TestPluginBase):
    package = 'q2_health_index.tests'

    # Minimal cases (only one H/N species present in the abundance table)

    def test_gmhi_predict_minimal_both_sp(self):
        table_file = self.get_data_path("input/abundances"
                                        "/minimal_data_both_sp.qza")
        table = qiime2.Artifact.load(table_file)
        res = health_index.actions.gmhi_predict(table=table)
        gmhi = pd.to_numeric(res[0].view(pd.Series))
        gmhi_exp = pd.read_csv(
            self.get_data_path("expected/minimal_data_both_sp.tsv"),
            sep='\t', index_col=0, header=0, squeeze=True)
        pdt.assert_series_equal(
            gmhi, gmhi_exp, check_dtype=False, check_index_type=False,
            check_series_type=False, check_names=False)

    def test_gmhi_predict_minimal_both_false(self):
        with self.assertRaisesRegex(AssertionError,
                                    "Could not find healthy species"):
            table_file = self.get_data_path("input/abundances"
                                            "/minimal_data_both_false.qza")
            table = qiime2.Artifact.load(table_file)
            health_index.actions.gmhi_predict(table=table)

    def test_gmhi_predict_minimal_one_h_sp(self):
        with self.assertRaisesRegex(AssertionError,
                                    "Could not find non-healthy species"):
            table_file = self.get_data_path("input/abundances"
                                            "/minimal_data_one_h_sp.qza")
            table = qiime2.Artifact.load(table_file)
            health_index.actions.gmhi_predict(table=table)

    def test_gmhi_predict_minimal_one_nh_sp(self):
        with self.assertRaisesRegex(AssertionError,
                                    "Could not find healthy species"):
            table_file = self.get_data_path("input/abundances"
                                            "/minimal_data_one_nh_sp.qza")
            table = qiime2.Artifact.load(table_file)
            health_index.actions.gmhi_predict(table=table)

    # Feature table with full taxonomy (real-world scenario)

    def test_gmhi_predict_full_taxonomy(self):
        table_file = self.get_data_path("input/abundances"
                                        "/full_taxonomy_mock_feature_table.qza")
        table = qiime2.Artifact.load(table_file)
        res = health_index.actions.gmhi_predict(table=table)
        gmhi = pd.to_numeric(res[0].view(pd.Series))
        gmhi_exp = pd.read_csv(
            self.get_data_path("expected/mock_data_only_species.tsv"),
            sep='\t', index_col=0, header=0, squeeze=True)
        pdt.assert_series_equal(
            gmhi, gmhi_exp, check_dtype=False, check_index_type=False,
            check_series_type=False, check_names=False)

    # Basic examples (dataset from Gupta et al. 2020)

    def test_gmhi_predict_4347_final(self):
        table_file = self.get_data_path("input/abundances"
                                        "/4347_final_relative_abundances.qza")
        table = qiime2.Artifact.load(table_file)
        res = health_index.actions.gmhi_predict(table=table)
        gmhi = pd.to_numeric(res[0].view(pd.Series))
        gmhi_exp = pd.read_csv(
            self.get_data_path("expected/4347_final_gmhi_Python.tsv"),
            sep='\t', index_col=0, header=0, squeeze=True)
        pdt.assert_series_equal(
            gmhi, gmhi_exp, check_dtype=False, check_index_type=False,
            check_series_type=False, check_names=False)

    def test_gmhi_predict_viz_4347_final(self):
        table_file = self.get_data_path("input/abundances"
                                        "/4347_final_relative_abundances.qza")
        metadata_file = self.get_data_path("input/metadata"
                                           "/4347_final_metadata.tsv")
        table = qiime2.Artifact.load(table_file)
        metadata = qiime2.Metadata.load(metadata_file)
        res = health_index.actions.gmhi_predict_viz(table=table,
                                                      metadata=metadata)
        gmhi = pd.to_numeric(res[0].view(pd.Series))
        gmhi_exp = pd.read_csv(
            self.get_data_path("expected/4347_final_gmhi_Python.tsv"),
            sep='\t', index_col=0, header=0, squeeze=True)
        pdt.assert_series_equal(
            gmhi, gmhi_exp, check_dtype=False, check_index_type=False,
            check_series_type=False, check_names=False)


if __name__ == '__main__':
    unittest.main()
