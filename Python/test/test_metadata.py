import unittest
import pandas as pd
import atlasapprox_disease as aad


class TestMetadata(unittest.TestCase):
    def setUp(self):
        """
        Set up the API instance before each test.
        """
        # os.environ["ATLASAPPROX_DISEASE_BASEURL"] = "http://127.0.0.1:5000/v1/"
        self.api = aad.API()  # Initialize the API instance

    def test_metadata_monocyte_lung(self):
        """
        Test metadata retrieval for monocyte-related lung diseases.
        """
        result = self.api.metadata(cell_type="monocyte", tissue="lung")
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0, "Metadata query returned empty DataFrame.")

        # Assert that all cell types contain 'monocyte'
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

        # Assert that all sex values are 'female'
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

        # Validate filtering correctness
        for _, row in result.iterrows():
            self.assertIn("kidney", row["disease"].lower())
            self.assertIn("epithelial", row["cell_type"].lower())
            self.assertEqual(row["sex"].lower(), "female")

