JavaScript
=========

The JavaScript interface can be used to access the atlasapprox-disease API from node.js or a web page.

Quick Start
-----------

.. code-block:: javascript

   let average_expression;
      (async () => {
      average_expression = await atlasapprox_disease.average(
         disease = "covid",
         tisse = "lung",
      );
      console.log(average_expression);
   })();

Installation
------------

You can use npm to install the atlasapprox-disease package:

.. code-block:: bash

   npm install atlasapprox-disease

Getting Started
---------------

Import or require the API object:

.. code-block:: javascript

   api = require('atlasapprox-disease');

Reference API
-------------

The complete API reference for the JavaScript interface is outlined below.

metadata
++++++++

**Description**:

Retrieves metadata records from the atlasapprox-disease API. Each record represents a unique combination of dataset, cell type, tissue, disease condition, sex, and developmental stage that meets the query criteria.

**Parameters**:

- ``disease`` *(optional)* – Filter by disease name (e.g., "COVID").
- ``cell_type`` *(optional)* – Filter by cell type (e.g., "T cell").
- ``tissue`` *(optional)* – Filter by tissue (e.g., "lung").
- ``sex`` *(optional)* – Filter by sex (e.g., "male", "female").
- ``development_stage`` *(optional)* – Filter by developmental stage (e.g., "adult").

**Returns**:

A promise that resolves to an object containing the metadata records.

differential_cell_type_abundance
++++++++++++++++++++++++++++++++

**Description**:

Retrieves differential cell type abundance across conditions such as disease, tissue, sex, or developmental stage. This method enables users to compare cell type proportions between selected conditions.

**Parameters**:

- ``differential_axis`` *(default: "disease")* – The axis for comparison (e.g., "disease", "sex").
- ``disease`` *(optional)* – Filter by disease name (e.g., "covid").
- ``cell_type`` *(optional)* – Filter by cell type (e.g., "macrophage").
- ``tissue`` *(optional)* – Filter by tissue (e.g., "lung").
- ``sex`` *(optional)* – Filter by sex (e.g., "male", "female").
- ``development_stage`` *(optional)* – Filter by developmental stage (e.g., "adult").

**Returns**:

A promise that resolves to an object containing the differential cell type abundance data.

differential_gene_expression
++++++++++++++++++++++++++++

**Description**:

Retrieves differentially expressed genes between a baseline condition and a specified state (e.g., disease vs. normal). By default, it identifies the top 10 up and down-regulated genes in each cell type across all datasets that match the filter criteria.

**Parameters**:

- ``differential_axis`` *(default: "disease")* – The axis for comparison (e.g., "disease", "sex").
- ``disease`` *(optional)* – Filter by disease name (e.g., "covid").
- ``cell_type`` *(optional)* – Filter by cell type (e.g., "macrophage").
- ``tissue`` *(optional)* – Filter by tissue (e.g., "lung").
- ``sex`` *(optional)* – Filter by sex (e.g., "male", "female").
- ``development_stage`` *(optional)* – Filter by developmental stage (e.g., "adult").
- ``top_n`` *(optional)* – Number of top differentially expressed genes to return (default: 10). Cannot be used with ``feature``.
- ``feature`` *(optional)* – The gene to query. Cannot be used with ``top_n``.
- ``method`` *(default: "delta_fraction")* – Method to calculate differential expression ("delta_fraction" or "ratio_average").

**Returns**:

A promise that resolves to an object containing the differential gene expression data.

highest_measurement
+++++++++++++++++++

**Description**:

Retrieves the top N cell types and tissue combinations with the highest expression of a given feature (gene) across multiple datasets. This helps identify the most highly expressing cell types for a gene of interest in different diseases and tissues.

**Parameters**:

- ``feature`` *(required)* – The gene to query.
- ``number`` *(optional)* – Number of highest expressing cell types to return (default: 10).

**Returns**:

A promise that resolves to an object containing the highest measurement data, ordered by expression level.

average
+++++++

**Description**:

Retrieves the average expression levels of one or more selected features (e.g., genes) across cell types, tissues, and diseases.

**Parameters**:

- ``features`` *(required)* – A comma-separated string or array of features (genes) to query.
- ``disease`` *(optional)* – Filter by disease (e.g., "covid").
- ``cell_type`` *(optional)* – Filter by cell type (e.g., "T cell").
- ``tissue`` *(optional)* – Filter by tissue (e.g., "lung").
- ``sex`` *(optional)* – Filter by sex (e.g., "male", "female").
- ``development_stage`` *(optional)* – Filter by developmental stage (e.g., "adult").
- ``unique_ids`` *(optional)* – The unique_ids from metadata results.
- ``include_normal`` *(optional)* – Include the corresponding normal condition if true (default: false).

**Returns**:

A promise that resolves to an object containing the average expression data.

fraction_detected
+++++++++++++++++

**Description**:

Retrieves the fraction of cells in which a given gene is detected across different cell types, tissues, and diseases. This provides an estimation of how commonly a gene is expressed in a given cell population.

**Parameters**:

- ``features`` *(required)* – A comma-separated string or array of features (genes) to query.
- ``disease`` *(optional)* – Filter by disease (e.g., "covid").
- ``cell_type`` *(optional)* – Filter by cell type (e.g., "T cell").
- ``tissue`` *(optional)* – Filter by tissue (e.g., "lung").
- ``sex`` *(optional)* – Filter by sex (e.g., "male", "female").
- ``development_stage`` *(optional)* – Filter by developmental stage (e.g., "adult").
- ``unique_ids`` *(optional)* – The unique_ids from metadata results.
- ``include_normal`` *(optional)* – Include the corresponding normal condition if true (default: false).

**Returns**:

A promise that resolves to an object containing the fraction detected data.

dotplot
+++++++

**Description**:

Retrieves both the average expression and fraction detected for a list of genes across different cell types, tissues, and diseases. This method is used for visualizing gene expression in a dot plot format, where dot size represents fraction detected and color represents average expression.

**Parameters**:

- ``features`` *(required)* – A comma-separated string or array of features (genes) to query.
- ``disease`` *(optional)* – Filter by disease (e.g., "covid").
- ``cell_type`` *(optional)* – Filter by cell type (e.g., "T cell").
- ``tissue`` *(optional)* – Filter by tissue (e.g., "lung").
- ``sex`` *(optional)* – Filter by sex (e.g., "male", "female").
- ``development_stage`` *(optional)* – Filter by developmental stage (e.g., "adult").
- ``unique_ids`` *(optional)* – The unique_ids from metadata results.
- ``include_normal`` *(optional)* – Include the corresponding normal condition if true (default: false).

**Returns**:

A promise that resolves to an object containing the dot plot data.