import unittest

import pandas as pd
import atlasapprox_disease as aad


class TestDifferentialCellTypeAbundance(unittest.TestCase):
    def setUp(self):
        """
        Set up the API instance before each test.
        """
        self.api = aad.API()  # Initialize the API instance
        
    def test_basic(self):
        result = self.api.differential_cell_type_abundance(
            tissue="lung",
            cell_type="T cell",
        )
        self.assertEqual(type(result), pd.DataFrame)
        self.assertGreater(len(result), 0)

    def test_default_all_results(self):
        """
        Test that no filters return all results.
        """
        result = self.api.differential_cell_type_abundance()
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0)  # Ensure the result is not empty
