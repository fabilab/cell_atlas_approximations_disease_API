# import unittest
# import pandas as pd
# import atlasapprox_disease as aad

# class TestHighestMeasurement(unittest.TestCase):
#     def setUp(self):
#         """
#         Set up the API instance before each test.
#         """
#         self.api = aad.API()

#     def test_highest_measurement(self):
#         """
#         Test that the function will return 10 rows of data by default
#         """
#         result = self.api.highest_measurement(
#             feature="KRAS",
#         )
        
#         self.assertIsInstance(result, pd.DataFrame)
#         self.assertTrue(
#             result.shape[0] == 10,
#             "By default, the return dataframe should have 10 rows."
#         )

    
    
        
        