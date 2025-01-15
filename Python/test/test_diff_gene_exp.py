import unittest
import atlasapprox_disease as aad
import pandas as pd


class TestDifferentialGeneExpression(unittest.TestCase):
    def setUp(self):
        """
        Set up the API instance before each test.
        """
        self.api = aad.API()  # Initialize the API instance

    def test_basic(self):
        result = self.api.differential_gene_expression(
            disease="COVID", cell_type="T cell", tissue="lung"
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0)
        self.assertIn("gene", result.columns)
        self.assertIn("metric", result.columns)

    def test_specific_feature(self):
        """
        Test fetching differential expression for specific features.
        """
        result = self.api.differential_gene_expression(
            disease="COVID-19", cell_type="T cell", tissue="lung", feature="IL6"
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0)
        self.assertIn("gene", result.columns)
        self.assertIn("metric", result.columns)
        self.assertTrue(all(result["gene"] == "IL6"))

