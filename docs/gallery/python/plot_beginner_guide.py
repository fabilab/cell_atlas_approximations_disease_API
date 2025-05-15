"""
.. _beginner-guide:

Beginner guide
==============

The `atlasapprox-disease <https://github.com/fabilab/cell_atlas_approximations_disease_API>`_ Python API provides access to over 600 disease-related single-cell datasets.
Currently, it includes datasets from the CELLxGENE Census as its initial source, covering diseases such as COVID-19, diabetes, acute kidney failure, and gastritis, along with metadata like cell type, developmental stage, and sex.
This API enables users to quickly explore cellular and gene expression patterns in disease contexts.

Follow this tutorial to get started with the basics of using the API.
"""

# %%
# Installation
# ------------
#
# (Optional) To ensure consistent dependencies, we recommend setting up a virtual environment:
#
# .. code-block:: bash
#
#     python -m venv ./venv
#     source ./venv/bin/activate
#
# Then, install the ``atlasapprox-disease`` package using ``pip``:
#
# .. code-block:: bash
#
#     pip install atlasapprox-disease

# %%
# Python quick start
# ------------------
#
# Below are 2 examples of common operations you can do with the atlasapprox_disease Python API:

# Import the package and initialise the API
import atlasapprox_disease

api = atlasapprox_disease.API()

# %%
# Querying cell metadata
# ^^^^^^^^^^^^^^^^^^^^^^
#
# The ``metadata`` function lets you explore cell metadata across datasets by applying filters
# on attributes like tissue, disease, or developmental stage.
# For example, the following filters cells from lung tissue at the adult stage:

api.metadata(
    tissue="lung",
    development_stage="adult"
)

# %%
# The output is a ``pandas.DataFrame`` with over 1400 unique combinations of cell types, diseases,
# and other columns such as ``sex``, ``cell_count`` and the ``dataset`` it comes from.
# You can get a quick overview of what cell types and conditions are available
# and use this information later for querying other API functions.

# %%
# Querying average gene expression
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The ``average`` function retrieves average gene expression across cell types, tissues, and diseases.
# For example, to query immune-related genes in COVID-19:

api.average(
    features="ACE2,TLR4,NLRP3,MBL2,IL6",
    disease="COVID-19"
)

# %%
# The output is a ``pandas.DataFrame`` with columns such as ``cell_type``, ``tissue_general``,
# ``disease``, ``dataset_id``, and the expression levels of the queried genes (in counts per 10k).
# This helps you explore gene activity in specific conditions and identify key genes
# for further analysis.

# %%
# Next steps
# ----------
#
# This tutorial introduced the basics of the ``atlasapprox-disease`` API.
# To learn more, explore additional functions like ``dotplot`` for visualizing gene expression,
# or query ``differential gene expression`` data.
#
# Visit the `official documentation <https://cell-atlas-approximations-disease-api.readthedocs.io/en/latest/python/index.html>`_
# for further details.
