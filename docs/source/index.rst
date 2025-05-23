.. atlasapprox-disease documentation master file, created by
   sphinx-quickstart on Fri Jan 17 14:05:38 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

API for disease cell atlas approximations
===============================================

Disease cell atlases are single-cell omics datasets that provide insights into diseased tissues across multiple organs and conditions. A disease `cell atlas approximation <https://chanzuckerberg.com/science/programs-resources/single-cell-biology/data-insights/light-and-scalable-statistical-approximations-of-cell-atlases>`_ is a lightweight, lossy compression of such a cell atlas, retaining key disease-specific information while reducing data size and complexity. This API, initially sourced from the CellxGene Census, enables researchers to explore cell type and gene expression data between disease states, sexes, and other conditions. It helps answer biological questions such as:

- *What cell types are present in blood samples from adult patients with influenza across all datasets?*
- *How does the proportion of alveolar type 2 cells in the lung change between healthy individuals and those with COVID-19*
- *What are the top 20 most differentially expressed genes in kidney-related disease?*

Version
-------
The most recent version of the API is **v1**.

Interfaces
----------
There are multiple ways to access atlas approximations programmatically:

.. toctree::
   :maxdepth: 1

   Python <python/index>
   REST (language-agnostic) <rest/index>
   JavaScript <js/index>

Authors
-------
This project is developed and maintained by:
**Fabio Zanini**, **Ying Xu**

Affiliated with **[Fabilab](https://fabilab.org)**.

Citation
--------

- **Xu et al.** (2024). Lightweight and scalable approximations democratise access to single cell atlases. *bioRxiv*. `doi:10.1101/2024.01.03.573994 <https://www.biorxiv.org/content/10.1101/2024.01.03.573994v1>`_.

Data sources
-----------

CZ CELLxGENE: Discover (Census). CZ CELLxGENE Discover: A single-cell data platform for scalable exploration, analysis and modeling of aggregated data CZI Single-Cell Biology, et al. bioRxiv 2023.10.30; doi: https://doi.org/10.1101/2023.10.30.563174