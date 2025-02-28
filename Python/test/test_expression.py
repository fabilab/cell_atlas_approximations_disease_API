import unittest
import pandas as pd
import atlasapprox_disease as aad

# Unit tests for measurement-related functions: average, fraction_detected, and dotplot.
class TestMeasurement(unittest.TestCase):
    def setUp(self):
        """
        Set up the API instance before each test.
        """
        self.api = aad.API()

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
        self.assertTrue(
            result["cell_type"].str.contains("endothelial", case=False, na=False).all(),
            "Not all cell types contain 'endothelial' in the result."
        )

    def test_average_filtered_by_disease_with_normal(self):
        """
        Test that when `include_normal=True`, the results includes both normal and disease data.
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
        self.assertTrue(
            result["disease"].str.contains("normal|covid", case=False, na=False).all(),
            "Results do not contain both 'normal' and 'covid'."
        )

    def test_fraction_sex_column_presence(self):
        """
        Verify the presence of the 'sex' column when filtering by sex.
        """
        result_with_sex = self.api.fraction_detected(
            features="INS,GCK,MAFA,PECAM1",
            disease="Diabetes",
            cell_type="endothelial",
            sex="male",
            development_stage="adult",
        )

        self.assertIn("sex", result_with_sex.columns, "Expected 'sex' column when filtering by sex.")
        
        result_without_sex = self.api.fraction_detected(
            features="INS,GCK,MAFA,PECAM1",
            disease="Diabetes",
            cell_type="endothelial",
            development_stage="adult",
        )

        self.assertNotIn("sex", result_without_sex.columns, "Unexpected 'sex' column when no sex filter is applied.")

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

    def test_dotplot_with_unique_ids(self):
        """
        Test that filtering by unique_id in dotplot returns expected values.
        """
        expected_genes = {"APOL1", "MYH9", "HNF1B"}
        expected_fraction_genes = {f"fraction_{gene}" for gene in expected_genes}
        expected_values = {
            "dataset_id": "0b75c598-0893-4216-afe8-5414cab7739d",
            "cell_type": "B cell",
            "tissue_general": "kidney",
            "disease": "acute kidney failure",
            "development_stage_general": "adult",
            "sex": "female",
        }
        result = self.api.dotplot(
            features="APOL1,MYH9,HNF1B",
            unique_ids="c8d24ef26af50f3c860ec433584bd9ad"
        )

        self.assertIsInstance(result, pd.DataFrame)
        assert expected_genes.issubset(result.columns), \
            f"Result is missing expected gene columns: {expected_genes}"

        assert expected_fraction_genes.issubset(result.columns), \
            f"Result is missing expected fraction-detected columns: {expected_fraction_genes}"

        for key, expected in expected_values.items():
            self.assertIn(key, result.columns, f"Missing column: {key}")
            self.assertEqual(result.iloc[0][key], expected, f"Mismatch in {key}: expected {expected}, got {result.iloc[0][key]}")
