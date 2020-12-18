import unittest
from warnings import filterwarnings

import biom
import pandas as pd
import pandas.util.testing as pdt
import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import health_index

from q2_health_index._utilities import (_load_file,
                                        _load_metadata,
                                        _load_and_validate_species,
                                        _validate_metadata_is_superset,
                                        _validate_and_extract_healthy_states,
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
        with self.assertRaisesRegex(FileNotFoundError, "No such file or directory"):
            infile = self.get_data_path("input/species/do-not-exists.txt")
            _load_and_validate_species(healthy_species_fp=infile)

    def test_species_load_non_healthy_wrong(self):
        with self.assertRaisesRegex(FileNotFoundError, "No such file or directory"):
            infile = self.get_data_path("input/species/do-not-exists.txt")
            _load_and_validate_species(non_healthy_species_fp=infile)

    def test_species_load_healthy_empty(self):
        with self.assertRaisesRegex(ValueError, "Healthy species list is empty!"):
            infile = self.get_data_path("input/species/empty_species.txt")
            _load_and_validate_species(healthy_species_fp=infile)

    def test_species_load_non_healthy_empty(self):
        with self.assertRaisesRegex(ValueError, "Non-healthy species list is empty!"):
            infile = self.get_data_path("input/species/empty_species.txt")
            _load_and_validate_species(non_healthy_species_fp=infile)

    def test_species_load_healthy_fake(self):
        infile = self.get_data_path("input/species/fake_MH_species.txt")
        healthy, non_healthy = _load_and_validate_species(healthy_species_fp=infile)
        self.assertListEqual(healthy, ['s__fake_1', 's__fake_2'])
        self.assertListEqual(non_healthy, NON_HEALTHY_SPECIES_DEFAULT)

    def test_species_load_non_healthy_fake(self):
        infile = self.get_data_path("input/species/fake_MN_species.txt")
        healthy, non_healthy = _load_and_validate_species(non_healthy_species_fp=infile)
        self.assertListEqual(non_healthy, ['s__fake_non_1', 's__fake_non_2', 's__fake_non_3'])
        self.assertListEqual(healthy, HEALTHY_SPECIES_DEFAULT)

    def test_species_load_fake(self):
        infile1 = self.get_data_path("input/species/fake_MH_species.txt")
        infile2 = self.get_data_path("input/species/fake_MN_species.txt")
        healthy, non_healthy = _load_and_validate_species(infile1, infile2)
        self.assertListEqual(healthy, ['s__fake_1', 's__fake_2'])
        self.assertListEqual(non_healthy, ['s__fake_non_1', 's__fake_non_2', 's__fake_non_3'])

    def test_metadata_not_provided(self):
        with self.assertRaisesRegex(ValueError, "Metadata parameter not provided!"):
            _load_metadata(None)

    def test_metadata_load_simple(self):
        infile = self.get_data_path("input/metadata/simple_metadata.tsv")
        metadata = _load_metadata(qiime2.Metadata.load(infile))
        metadata_exp = pd.read_csv(infile, sep='\t')
        self.assertListEqual(list(metadata.columns), ['Age', 'Healthy'])
        self.assertListEqual(list(metadata.index), list(metadata_exp['sample-id']))
        self.assertListEqual(list(metadata.Healthy), list(metadata_exp['Healthy']))
        self.assertListEqual(list(metadata.Age), list(metadata_exp['Age']))

    def test_metadata_validate_simple(self):
        table_file = self.get_data_path("input/abundances/simple_relative_abundances.qza")
        metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
        metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
        table = qiime2.Artifact.load(table_file)
        _validate_metadata_is_superset(metadata, table.view(biom.Table))

    def test_metadata_validate_simple_wrong(self):
        with self.assertRaisesRegex(ValueError, "Missing samples in metadata: {'MOCK-"):
            table_file = self.get_data_path("input/abundances/dada2_table.qza")
            metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
            metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
            table = qiime2.Artifact.load(table_file)
            _validate_metadata_is_superset(metadata, table.view(biom.Table))

    def test_healthy_states_correct_simple(self):
        metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
        metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
        hs, ns = _validate_and_extract_healthy_states(metadata, 'Healthy', 'Y', 'N,Sick')
        self.assertListEqual(hs, ['Y'])
        self.assertListEqual(ns, ['N', 'Sick'])
        _validate_and_extract_healthy_states(metadata, 'Healthy', 'Y', 'rest')
        _validate_and_extract_healthy_states(metadata, 'Healthy', 'rest', 'N,Sick')

    def test_healthy_states_parameters_not_provided(self):
        metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
        metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
        with self.assertRaisesRegex(ValueError, "healthy_column parameter not provided!"):
            _validate_and_extract_healthy_states(metadata, None, 'Y', 'N,Sick')
        with self.assertRaisesRegex(ValueError, "healthy_states parameter not provided!"):
            _validate_and_extract_healthy_states(metadata, 'Col', None, 'Par')
        with self.assertRaisesRegex(ValueError, "non_healthy_states parameter not provided!"):
            _validate_and_extract_healthy_states(metadata, 'Col', 'Par', None)

    def test_healthy_states_wrong_healthy_column(self):
        with self.assertRaisesRegex(ValueError, "\'Col\' is not a column in your metadata."):
            metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
            metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
            _validate_and_extract_healthy_states(metadata, 'Col', 'Par', 'Par')

    def test_healthy_states_both_states_rest(self):
        with self.assertRaisesRegex(ValueError, f'healthy_states and non_healthy_states '
                                                f'parameters cannot be equal.'):
            metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
            metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
            _validate_and_extract_healthy_states(metadata, 'Healthy', 'rest', 'rest')

    def test_healthy_states_both_states_equal(self):
        with self.assertRaisesRegex(ValueError, f'healthy_states and non_healthy_states '
                                                f'parameters cannot be equal.'):
            metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
            metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
            _validate_and_extract_healthy_states(metadata, 'Healthy', 'N,Sick,N', 'Sick,N,Sick')

    def test_healthy_states_healthy_state_wrong(self):
        with self.assertRaisesRegex(ValueError, f'Healthy state \'Par1\' is not represented '
                                                f'by any members of \'Healthy\' column in metadata.'):
            metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
            metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
            _validate_and_extract_healthy_states(metadata, 'Healthy', 'Par1', 'Par2')

    def test_healthy_states_non_healthy_state_wrong(self):
        with self.assertRaisesRegex(ValueError, f'Non-healthy state \'Par2\' is not represented '
                                                f'by any members of \'Healthy\' column in metadata.'):
            metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
            metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
            _validate_and_extract_healthy_states(metadata, 'Healthy', 'Y', 'Par2')

    def test_healthy_states_non_healthy_state_wrong(self):
        with self.assertRaisesRegex(ValueError, f'Number of healthy and non-healthy state values '
                                                f'is not equal to the number of rows in metadata.'):
            metadata_file = self.get_data_path("input/metadata/simple_metadata.tsv")
            metadata = _load_metadata(qiime2.Metadata.load(metadata_file))
            _validate_and_extract_healthy_states(metadata, 'Healthy', 'Y', 'N')


class TestCalculateGmhi(TestPluginBase):
    package = 'q2_health_index.tests'

    def test_calculate_gmhi_4347_final(self):
        table_file = self.get_data_path("input/abundances/4347_final_relative_abundances.qza")
        metadata_file = self.get_data_path("input/metadata/4347_final_metadata.tsv")
        table = qiime2.Artifact.load(table_file)
        metadata = qiime2.Metadata.load(metadata_file)
        res = health_index.actions.calculate_gmhi(
            table=table,
            metadata=metadata,
            healthy_column='phenotype',
            healthy_states='Healthy',
            non_healthy_states='rest',
            healthy_species_fp=None,
            non_healthy_species_fp=None)
        gmhi = pd.to_numeric(res[0].view(pd.Series))
        gmhi_exp = pd.read_csv(
            self.get_data_path("expected/4347_final_gmhi.tsv"),
            sep='\t', index_col=0, header=0, squeeze=True)
        pdt.assert_series_equal(
            gmhi, gmhi_exp, check_dtype=False, check_index_type=False,
            check_series_type=False, check_names=False)


if __name__ == '__main__':
    unittest.main()
