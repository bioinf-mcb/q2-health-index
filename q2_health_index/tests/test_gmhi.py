import unittest
import pandas as pd
import pandas.util.testing as pdt

from warnings import filterwarnings

import qiime2

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import health_index
from q2_health_index._utilities import (_load_and_validate_species, _load_metadata,
                                        HEALTHY_SPECIES_DEFAULT, NON_HEALTHY_SPECIES_DEFAULT)

filterwarnings("ignore", category=UserWarning)
filterwarnings("ignore", category=RuntimeWarning)


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
            _load_and_validate_species(list_healthy=infile)

    def test_species_load_non_healthy_wrong(self):
        with self.assertRaisesRegex(FileNotFoundError, "No such file or directory"):
            infile = self.get_data_path("input/species/do-not-exists.txt")
            _load_and_validate_species(list_non_healthy=infile)

    def test_species_load_healthy_empty(self):
        with self.assertRaisesRegex(ValueError, "Healthy species list is empty!"):
            infile = self.get_data_path("input/species/empty_species.txt")
            _load_and_validate_species(list_healthy=infile)

    def test_species_load_non_healthy_empty(self):
        with self.assertRaisesRegex(ValueError, "Non-healthy species list is empty!"):
            infile = self.get_data_path("input/species/empty_species.txt")
            _load_and_validate_species(list_non_healthy=infile)

    def test_species_load_healthy_fake(self):
        infile = self.get_data_path("input/species/fake_MH_species.txt")
        healthy, non_healthy = _load_and_validate_species(list_healthy=infile)
        self.assertListEqual(healthy, ['s__fake_1', 's__fake_2'])
        self.assertListEqual(non_healthy, NON_HEALTHY_SPECIES_DEFAULT)

    def test_species_load_non_healthy_fake(self):
        infile = self.get_data_path("input/species/fake_MN_species.txt")
        healthy, non_healthy = _load_and_validate_species(list_non_healthy=infile)
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

    def test_metadata_simple(self):
        infile = self.get_data_path("input/metadata/simple_metadata.tsv")
        metadata = _load_metadata(qiime2.Metadata.load(infile))
        metadata_exp = pd.read_csv(infile, sep='\t')
        self.assertListEqual(list(metadata.columns), ['Age', 'Healthy'])
        self.assertListEqual(list(metadata.index), list(metadata_exp['sample-id']))
        self.assertListEqual(list(metadata.Healthy), list(metadata_exp['Healthy']))
        self.assertListEqual(list(metadata.Age), list(metadata_exp['Age']))

    def test_calculate_gmhi_4347_final(self):
        table_file = self.get_data_path("input/abundances/4347_final_relative_abundances.qza")
        metadata_file = self.get_data_path("input/metadata/4347_final_metadata.tsv")
        table = qiime2.Artifact.load(table_file)
        metadata = qiime2.Metadata.load(metadata_file)
        res = health_index.actions.calculate_gmhi(
                             table=table,
                             metadata=metadata,
                             healthy_species=None,
                             non_healthy_species=None)
        gmhi = pd.to_numeric(res[0].view(pd.Series))
        gmhi_exp = pd.read_csv(
            self.get_data_path("expected/4347_final_gmhi.tsv"),
            sep='\t', index_col=0, header=0, squeeze=True)
        pdt.assert_series_equal(
            gmhi, gmhi_exp, check_dtype=False, check_index_type=False,
            check_series_type=False, check_names=False)


class TestCalculateGmhi(TestPluginBase):
    package = 'q2_health_index.tests'


if __name__ == '__main__':
    unittest.main()
