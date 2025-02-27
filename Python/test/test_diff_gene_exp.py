import unittest
import pandas as pd
import atlasapprox_disease as aad

class TestDifferentialGeneExpression(unittest.TestCase):
    def setUp(self):
        self.api = aad.API() 

    def test_diff_gene_disease_fibroblast(self):
        """Test differential gene expression in fibroblasts (Diabetes vs. Normal)."""
        top_n = 5
        result = self.api.differential_gene_expression(
            differential_axis="disease",
            disease="diabetes",
            cell_type="fibroblast",
            top_n=top_n,
            method="delta_fraction",
        )

        self.assertIsInstance(result, pd.DataFrame, "Result is not a pandas DataFrame.")
        self.assertTrue(
            result["cell_type"].str.contains("fibroblast", case=False, na=False).all(),
            "Not all rows have 'fibroblast' in the cell_type column."
        )
    
    def test_diff_gene_by_sex(self):
        """Test differential gene expression between male and female (Liver Hepatocytes)."""
        result = self.api.differential_gene_expression(
            differential_axis="sex",
            tissue="liver",
            cell_type="hepatocyte",
            method="ratio_average",
        )

        self.assertIsInstance(result, pd.DataFrame, "Result is not a pandas DataFrame.")
        self.assertTrue(
            (result["differential_axis"] == "sex").all(),
            "Not all rows have 'sex' in the differential_axis column."
        )

        self.assertTrue(
            result["cell_type"].str.contains("hepatocyte", case=False, na=False).all(),
            "Not all rows have 'hepatocyte' in the cell_type column."
        )

    def test_diff_gene_disease_specific_gene(self):
        """Test differential gene expression for ACE2 in kidney epithelial cells (Diabetes vs. Normal)."""
        result = self.api.differential_gene_expression(
            differential_axis="disease",
            disease="Diabetes",
            cell_type="epithelial cell",
            tissue="kidney",
            feature="ACE2",
            method="ratio_average",
        )
        self.assertTrue(
            (result["gene"] == "ACE2").all(),
            "Not all rows have 'ACE2' as the gene."
        )

