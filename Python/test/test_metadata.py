import unittest
import pandas as pd
import atlasapprox_disease as aad


class TestMetadata(unittest.TestCase):
    def setUp(self):
        """
        Set up the API instance before each test.
        """
        self.api = aad.API()

    def test_metadata_monocyte_lung(self):
        """
        Test metadata retrieval for monocyte-related lung diseases.
        """
        result = self.api.metadata(cell_type="monocyte", tissue="lung")
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0, "Metadata query returned empty DataFrame.")

        self.assertTrue(
            result["cell_type"].str.contains("monocyte", case=False, na=False).all(),
            "Not all cell types contain 'monocyte' in the result."
        )

    def test_metadata_female_kidney(self):
        """
        Test metadata retrieval for female kidney-related conditions.
        """
        result = self.api.metadata(tissue="kidney", sex="female")
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0, "Metadata query returned empty DataFrame.")

        self.assertTrue(
            (result["sex"].str.lower() == "female").all(),
            "Not all results have 'female' as the sex."
        )

    def test_metadata_multiple_filters(self):
        """
        Test metadata retrieval with multiple filters: disease, cell type, and sex.
        """
        result = self.api.metadata(disease="kidney", cell_type="epithelial", sex="female")
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0, "Metadata query returned empty DataFrame.")

        for _, row in result.iterrows():
            self.assertIn("kidney", row["disease"].lower())
            self.assertIn("epithelial", row["cell_type"].lower())
            self.assertEqual(row["sex"].lower(), "female")

