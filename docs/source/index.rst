.. atlasapprox-disease documentation master file, created by
   sphinx-quickstart on Fri Jan 17 14:05:38 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

API for disease cell atlas approximations
===============================================

Cell Atlas Approximation (Disease) is a lightweight API designed for analyzing approximated single-cell transcriptomics data, using CellxGene Census as its initial data source.
The API enables researchers to address complex biological questions across multiple organs, cell types, and disease conditions such as:

- *What is the expression of a specific gene in a specific disease across all datasets?*
- *In COVID-19, what is the differential cell abundance of each cell type between normal and disease states?*
- *What are the top 20 most differentially expressed genes in kidney disease?*

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
**Xu et al.** (2024) Lightweight and scalable approximations democratise access to single cell atlases. `biorxiv <https://www.biorxiv.org/content/10.1101/2024.01.03.573994v1>`_.