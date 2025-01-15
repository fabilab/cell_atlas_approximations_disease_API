import os
import unittest
import pandas as pd

import atlasapprox_disease as aad


class TestBasic(unittest.TestCase):
    def setUp(self):
        """
        Set up the API instance before each test.
        """
        os.environ["ATLASAPPROX_DISEASE_BASEURL"] = "http://127.0.0.1:5000/v1/"
        self.api = aad.API()  # Initialize the API instance

    def test_metadata(self):
        result = self.api.metadata(disease="Flu")
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0)
        self.assertIn("cell_type", result.columns)
        
    def test_metadata_invalid_keyword(self):
        """
        Test metadata retrieval with an invalid tissue filter.
        """
        result = self.api.metadata(tissue="hello")
        self.assertIsInstance(result, pd.DataFrame)
        self.assertTrue(result.empty)

    def test_metadata_with_multiple_filters(self):
        result = self.api.metadata(
            disease="kidney",
            cell_type="epithelial",
            sex="female",
        )
        self.assertEqual(type(result), pd.DataFrame)
        for _, row in result.iterrows():
            self.assertIn("kidney", row["disease"].lower())
            self.assertIn("epithelial", row["cell_type"].lower()) 
            
    def test_highest_measurement(self):
        result = self.api.highest_measurement(
            feature = "IL6",
            number = 15
        )
        self.assertEqual(type(result), pd.DataFrame)
        self.assertEqual(len(result), 15)
        

    def test_average(self):
        result = self.api.average(
            features = "IL6,CCL2,ACE2",
            disease = "Covid",
            cell_type = "macrophage",
        )
        
        self.assertEqual(type(result), pd.DataFrame)
        self.assertGreater(len(result), 0)
        # Assert all values in the cell_type column contain the substring 'macrophage'
        self.assertTrue(result["cell_type"].str.contains("macrophage", case=False, na=False).all())
    
    def test_fraction_detected(self):
        result = self.api.fraction_detected(
            features = "IL6,CCL2,ACE2",
            disease = "Covid",
            cell_type = "macrophage",
        )
        
        self.assertEqual(type(result), pd.DataFrame)
        self.assertGreater(len(result), 0)
        # Assert all values in the cell_type column contain the substring 'macrophage'
        self.assertTrue(result["cell_type"].str.contains("macrophage", case=False, na=False).all())
        
    def test_dotplot(self):
        result = self.api.dotplot(
            features = "IL6,CCL2,ACE2",
            disease = "Covid",
            cell_type = "macrophage",
        )
        
        self.assertEqual(type(result), pd.DataFrame)
        self.assertGreater(len(result), 0)
        # Assert all values in the cell_type column contain the substring 'macrophage'
        self.assertTrue(result["cell_type"].str.contains("macrophage", case=False, na=False).all())
        
    