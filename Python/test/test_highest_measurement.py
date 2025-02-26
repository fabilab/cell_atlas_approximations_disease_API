import unittest
import pandas as pd
import atlasapprox_disease as aad

# This pytest is designed to test three different functions that are related to measurement
# average, fraction_detected, and dotplot
class TestHighestMeasurement(unittest.TestCase):
    def setUp(self):
        """
        Set up the API instance before each test.
        """
        self.api = aad.API()  # Initialize the API instance

    def test_highest_measurement(self):
        """
        Test that the function will return 10 rows of data by default
        """
        result = self.api.highest_measurement(
            feature="KRAS",
        )
        
        self.assertTrue(
            result.shape[0] == 10,
            "By default, the return dataframe should have 10 rows."
        )
    
    def test_highest_measurement_given_number(self):
        """
        Test that the function will return the exact number of row queried by user
        """
        result = self.api.highest_measurement(
            feature="KRAS",
            number=30,
        )
        
        self.assertTrue(
            result.shape[0] == 30,
            "By default, the return dataframe should have 30 rows."
        )

    
    
        
        