import os
import unittest
import pandas as pd
import atlasapprox_disease as aad  # Ensure your package is installed

class TestDifferentialCellTypeAbundance(unittest.TestCase):
    def setUp(self):
        self.api = aad.API()  # Initialize the API instance

    def test_diff_cell_covid_lung(self):
        """
        Test differential cell type abundance in lung between COVID-19 patients and healthy controls.
        """
        result = self.api.differential_cell_type_abundance(
            disease="covid",
            tissue="lung"
        )
        
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0, "No data returned for COVID lung analysis.")
        self.assertIn("delta_frac", result.columns)
        self.assertTrue(
            result["disease"].str.contains("covid", case=False, na=False).all(),
            "Not all rows contain 'covid' in the disease column."
        )
        self.assertTrue(
            result["tissue_general"].str.contains("lung", case=False, na=False).all(),
            "Not all rows contain 'lung' in the tissue_general column."
        )

    def test_diff_cell_multi_filters(self):
        """
        Test differential cell type abundance with multiple filters.
        """
        result = self.api.differential_cell_type_abundance(
            tissue="lung",
            cell_type="macrophage",
            sex="female"
        )

        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0, "No data returned for macrophages in female lung.")
        
        # Ensure the filter conditions are met
        self.assertTrue(
            all(result["cell_type"].str.contains("macrophage", case=False, na=False)),
            "Result contains cell types other than macrophages."
        )
        self.assertTrue(
            all(result["sex"] == "female"),
            "Result contains sex values other than 'female'."
        )

    def test_diff_tcells_by_sex(self):
        """
        Test differential cell type abundance of T cells in blood between males and females.
        """
        result = self.api.differential_cell_type_abundance(
            differential_axis="sex",
            cell_type="T cell",
            tissue="blood"
        )

        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0, "No data returned for T cell abundance by sex.")
        self.assertIn("sex", result.columns)