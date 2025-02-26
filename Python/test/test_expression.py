import unittest
import pandas as pd
import atlasapprox_disease as aad

# This pytest is designed to test three different functions that are related to measurement
# average, fraction_detected, and dotplot
class TestMeasurement(unittest.TestCase):
    def setUp(self):
        """
        Set up the API instance before each test.
        """
        self.api = aad.API()  # Initialize the API instance

    def test_average_filtered_by_cell_type(self):
        """
        Test that filtering by cell_type only returns relevant cell types.
        """
        result = self.api.average(
            features="INS,GCK,MAFA,PECAM1",
            disease="Diabetes",
            cell_type="endothelial",
            development_stage="adult",
        )

        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0, "Query returned empty DataFrame.")

        # Ensure only endothelial cells are returned
        self.assertTrue(
            result["cell_type"].str.contains("endothelial", case=False, na=False).all(),
            "Not all cell types contain 'endothelial' in the result."
        )

    def test_average_filtered_by_disease_with_normal(self):
        """
        Test that when `include_normal=True`, the results contain both normal and disease data.
        """
        result = self.api.average(
            features="CD19, CD68, COL1A1",
            disease="covid",
            cell_type="T cell",
            sex="male",
            development_stage="adult",
            include_normal=True,
        )

        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0, "Query returned empty DataFrame.")

        # Ensure both 'normal' and 'covid' appear in the disease column
        self.assertTrue(
            result["disease"].str.contains("normal|covid", case=False, na=False).all(),
            "Results do not contain both 'normal' and 'covid'."
        )

    def test_fraction_sex_column_presence(self):
        """
        Test that if `sex` is specified, the 'sex' column exists; otherwise, it does not.
        """
        # Query with sex filter
        result_with_sex = self.api.fraction_detected(
            features="INS,GCK,MAFA,PECAM1",
            disease="Diabetes",
            cell_type="endothelial",
            sex="male",
            development_stage="adult",
        )

        self.assertIn("sex", result_with_sex.columns, "Column 'sex' should be present when filtering by sex.")

        # Query without sex filter
        result_without_sex = self.api.fraction_detected(
            features="INS,GCK,MAFA,PECAM1",
            disease="Diabetes",
            cell_type="endothelial",
            development_stage="adult",
        )

        self.assertNotIn("sex", result_without_sex.columns, "Column 'sex' should NOT be present when no sex filter is used.")

    def test_dotplot_gene_columns_exist(self):
        """
        Test that the result contains the expected gene columns in both average expression and fraction.
        """
        expected_genes = {"INS", "GCK", "MAFA", "PECAM1"}
        expected_genes_fraction = {f"fraction_{gene}" for gene in expected_genes}

        result = self.api.dotplot(
            features="INS,GCK,MAFA,PECAM1",
            disease="Diabetes",
            cell_type="endothelial",
            development_stage="adult",
        )

        self.assertTrue(
            expected_genes.issubset(result.columns),
            f"Result is missing expected gene columns:{expected_genes}"
        )
        
        self.assertTrue(
            expected_genes_fraction.issubset(result.columns),
            f"Result is missing expected fraction-detected columns:{expected_genes_fraction}"
        )