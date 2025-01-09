import pytest
import unittest

import atlasapprox_disease as aad


class TestBasic(unittest.TestCase):
    def test_init(self):
        api = aad.API()

    def test_metadata(self):
        api = aad.API()
        result = api.metadata(disease="Crohn")
        print("reulst type is---------")
        print(type(result))
        self.assertEqual(type(result), list)
        self.assertTrue(len(result) > 0)
        
    def test_metadata_invalid_keyword(self):
        api = aad.API()
        result = api.metadata(tissue="hello")
        
        self.assertEqual(type(result), list)
        self.assertTrue(len(result) == 0)
    
    def test_metadata_with_multiple_filters(self):
        api = aad.API()
        result = api.metadata(
            disease="kidney",
            cell_type="epithelial",
            sex="female",
        )
        # verify the response type
        self.assertEqual(type(result), list)
        # check that the response is not empty
        self.assertTrue(len(result) > 0)
        for record in result:
            assert "disease" in record
            assert "kidney" in record["disease"].lower()