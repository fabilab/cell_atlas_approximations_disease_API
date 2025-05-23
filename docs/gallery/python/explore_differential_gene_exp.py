"""
.. _differential-gene-exp:

Differential gene expression analysis
=====================================

This tutorial showcases one of the basic ways to perform differential gene expression analysis using the atlasapprox-disease API.
You will use the ``metadata`` and ``differential_gene_expression`` functions to identify datasets for a specific cell type, analyze gene expression changes in a disease context, and identify frequently occurring differentially expressed genes across datasets.
The tutorial uses memory B cells as an example, but you can apply the code to any cell type, disease, or tissue of interest using the API's many features.

# sphinx_gallery_thumbnail_path = '_static/differential_gene_exp.png'

"""


# %%
# Contents
# --------
#
# - Overview metadata with filters
# - Perform differential gene expression analysis for a specific cell type across datasets
# - Find frequently occurring differentially expressed genes
# - Tips for further exploration

# %%
# Installation
# ------------
#
# Install the required packages using `pip`:
#
# .. code-block:: bash
#
#     pip install atlasapprox-disease pandas

# %%
# Import libraries and initialize the API
# ---------------------------------------
#
# Import the necessary libraries
import atlasapprox_disease as aad
import pandas as pd

# Initialize the API
api = aad.API()

# %%
# Overview datasets with cell type-specific data
# ---------------------------------------------
#
# One way to start is to use the ``metadata`` function to get an overview of all the data relevant to what you want to explore, such as a cell type, disease, or tissue.
# In this example, we will focus on datasets related to memory B cells as a simple starting point:

cell_metadata = api.metadata(cell_type="memory B cell")

# Display the result
cell_metadata

#%%
# The DataFrame contains 186 rows, each representing a unique combination of metadata attributes (e.g. tissue, disease, sex, and development stage) involving memory B cells.


# %%
# To see the full list of unique diseases without truncation:
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

cell_metadata.disease.unique()

#%%
# As shown, there is a variety of diseases involving memory B cell data, e.g., COVID-19, post-COVID-19 disorder, breast carcinoma, and Crohn disease, which you can explore further.
# For example, you can select a disease like COVID-19 to perform differential gene expression analysis on memory B cells, as demonstrated in the following sections:

# %%
# Perform differential gene expression analysis for memory B cells in COVID-19
# ----------------------------------------------------------------------------
#
# To understand how memory B cells respond to COVID-19, query the top 10 up- and down-regulated genes (20 in total) across all datasets with diseased and normal conditions.
# This analysis identifies genes with the most significant expression changes in COVID-19 compared to healthy samples.

df_genes = api.differential_gene_expression(
    differential_axis = "disease",
    disease="covid",
    cell_type="memory B cell",
    top_n=10  # Top 10 up and down-regulated genes to query
)

# Display the results
df_genes

# %%
# The resulting DataFrame lists the top 10 up- and down-regulated genes for memory B cells in COVID-19 across all relevant datasets.
# Key columns include gene, regulation, expression and metric (fold change).
# Up-regulated genes may indicate activation of immune memory or antibody production pathways in response to COVID-19, while down-regulated genes could suggest suppression of other functions.
# Since the query includes multiple datasets and tissues, variations in gene expression may reflect dataset-specific or tissue-specific differences.


# %%
# Find frequently occurring differentially expressed genes
# --------------------------------------------------------
#
# Since memory B cells are present in multiple datasets, identify which genes appear most frequently as top differentially expressed genes across these datasets.
# This analysis highlights genes consistently affected by COVID-19 in memory B cells.

# Count the frequency of up-regulated genes across datasets
up_gene_counts = df_genes[df_genes["regulation"] == "up"]["gene"].value_counts()

# Display the results
print("Frequency of up-regulated genes across datasets:")
print(up_gene_counts)

# %%
# The output shows the frequency of up-regulated genes across datasets, for example, XAF1 appearing 4 times, NFKBID 3 times, MX1 3 times, and so on.
# This is how you can use the API to identify genes that frequently appear as top differentially expressed genes in your analysis.
# You can also explore down-regulated genes or analyze other diseases to compare results across different conditions.

# %%
# Examples for further exploration
# --------------------------------
#
# 1. Analyze down-regulated genes: Repeat the frequency analysis for down-regulated genes to identify consistently suppressed pathways.
#
#    .. code-block:: python
#
#        down_gene_counts = df_genes[df_genes["regulation"] == "down"]["gene"].value_counts()
#        print(down_gene_counts)
#
# 2. Explore other diseases: Use the diseases from the metadata (e.g., influenza) to compare memory B cell responses across conditions.
#
#    .. code-block:: python
#
#        df_influenza = api.differential_gene_expression(
#            disease="influenza",
#            cell_type="memory B cell",
#            top_n=10
#        )
#
# 3. Explore specific tissues: Query differential gene expression for a specific tissue (e.g., kidney) to analyze expression changes across all diseases and cell types in that tissue.
#
#    .. code-block:: python
#
#        df_kidney = api.differential_gene_expression(
#            tissue="kidney",
#            top_n=10
#        )
#        print(df_kidney)

# %%
# Next steps
# ----------
#
# This tutorial introduced differential gene expression analysis with the atlasapprox-disease API.
# To learn more, explore additional functions like average to retrieve gene expression levels, or dotplot for visualizing expression patterns.
#
# Visit the official documentation <https://cell-atlas-approximations-disease-api.readthedocs.io/en/latest/python/index.html>
# for further details.
