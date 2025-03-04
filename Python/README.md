Python interface to disease-related cell atlas approximations
=============================================

`atlasapprox_disease` is a Python interface designed to provide easy access to disease-related cell atlas data. Built using the [cellxgene census](https://github.com/chanzuckerberg/cellxgene-census) as the initial data source and powered by the [scquill approximation algorithm](https://github.com/fabilab/scquill/), this package enables researchers to explore thousands of datasets efficiently. 

The package enables researchers to address complex biological questions across multiple organs, cell types, and disease conditions such as:

- *What is the expression of a specific gene in a specific disease across all datasets?*
- *In COVID-19, what is the differential cell abundance of each cell type between normal and disease states?*
- *What are the top 20 most differentially expressed genes in kidney disease?*

---

## Features

- **Differential analysis**:
    - Compare differential gene expression across diseases, tissues, sexes, or conditions.
    - Retrieve differential cell type abundance between disease states.
- **Expression & Abundance Queries:**:
    - Retrieve average gene expression levels across multiple datasets
    - Compute the fraction of cells expressing specific genes.
    - Identify the highest-expressing cell types for a given gene across diseases and datasets.
- **Metadata Access**:
    - Retrieve dataset metadata, including cell types, tissues, diseases, sex, development stage.
- **Visualization:**:
    - Generate dot plots to visualize gene expression and detection frequencies.

## Installation

**To install the package, run:**
```bash
pip install atlasapprox-disease
```

## Quick start
**Initialise the API**
```python
import atlasapprox_disease as aad
api = aad.API()
```

**Retrieve Metadata**
```python
metadata = api.metadata()
print(metadata.head())
```

**Query differential gene expression**
```python
diff_exp = api.differential_gene_expression(
    differential_axis="disease",
    disease="Diabetes",
    top_n=10,  # Get top 10 upregulated & downregulated genes
    method="delta_fraction"
)
print(diff_exp)
```

**Generate average gene expression**
```python
average_exp = api.average(
    features="INS,GCK,MAFA,PECAM1",
    disease="Diabetes",
    include_normal=True
)
print(average_exp)
```

For more detailed tutorials and examples, please visit the [API github repo](https://github.com/fabilab/cell_atlas_approximations_disease_API).