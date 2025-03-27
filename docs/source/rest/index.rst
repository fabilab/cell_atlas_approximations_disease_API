REST
====
Cell atlas approximations (disease) are designed to be machine-readable and programming language-independent. This is achieved via a RESTful API.

The current version of the RESTful API is **v1**.

Quick Start
-----------
.. tabs::

   .. tab:: **Python**

      .. code-block:: python

         import requests
         import pandas as pd
         response = requests.get(
             'https://api-disease.atlasapprox.org/v1/metadata',
             params=dict(disease='COVID'),
         )
         print(pd.DataFrame(response.json()))

Getting started
---------------
- The API accepts both **GET** and **POST** requests.
- The API returns **JSON** data.
- No aliases for names (e.g.disease, genes) are supported yet: please double check your spelling.
- If you can use one of the languge-dedicated APIs (e.g. the Python API), please do so instead of using the REST API. Language-specific packages use caching to reduce load on our servers and also give you faster answers, so it’s a win-win.

Reference API
-------------
The complete API reference is outlined below. Each endpoint is described with its parameters and expected responses.  
For instance, to query the average expression of a gene list, the endpoint is ``/average`` and the full URL would be:  
``https://api-disease.atlasapprox.org/v1/average``


Metadata
++++++++
**Endpoint**: ``/metadata``

**Method**: ``GET``  

**Description**:

Retrieve metadata records from the source database that match the specified filters. Each returned record represents a unique combination of dataset, cell type, tissue, disease condition, sex, and developmental stage that meets the query criteria.

**Parameters**:

- ``disease`` *(optional)* – Filter by disease.  
- ``cell_type`` *(optional)* – Filter by cell type.  
- ``tissue`` *(optional)* – Filter by tissue.  
- ``sex`` *(optional)* – Filter by sex.  
- ``development_stage`` *(optional)* – Filter by developmental stage.

**Returns: A dict with the following key-value pairs:**

- ``unique_ids`` – The unique id that is specific for this combination of returned data.  
- ``dataset_id`` – Dataset ID of the where the cells come from.
- ``cell_type`` – The type of the cells as annotated in the dataset.  
- ``tissue_general`` – The general tissue category where the cells were observed.  
- ``disease`` – The disease condition associated with the dataset.  
- ``sex`` – The biological sex of the sample donor (e.g., male, female).  
- ``development_stage_general`` – The general developmental stage (e.g., adult, child).  
- ``cell_count`` – The number of recorded cells in the dataset that match the filter criteria.  

.. note::

   - Each object in the returned list represents a unique **metadata combination** that satisfies the applied filters.
   - The results are **split by cell type and sex**, meaning that if multiple cell types or sexes exist within a dataset, each combination is listed as a separate entry.
   - The same ``dataset_id`` may appear multiple times, as it can contain multiple cell types or sex-based subgroupings.

Differential cell type abundance
++++++++++++++++++++++++++++++++
**Endpoint**: ``/differential_cell_type_abundance``

**Method**: ``POST``  

**Description**:

Retrieve differential cell type abundance across conditions such as disease, tissue, sex, or developmental stage. This endpoint enables users to compare cell type proportions between selected conditions.

**Parameters**:

- ``differential_axis`` *(default: disease)* – The axis along which differential abundance is calculated, i.e. **disease** = disease condition vs normal condition, **sex** = male vs female, 
- ``disease`` *(optional)* – Filter by disease condition.  
- ``cell_type`` *(optional)* – Filter by cell type.  
- ``tissue`` *(optional)* – Filter by tissue.  
- ``sex`` *(optional)* – Filter by sex.  
- ``development_stage`` *(optional)* – Filter by developmental stage.  

**Returns: A dict with the following key-value pairs:**  

- ``dataset_id`` – Dataset ID of the cells that satisfies the filter conditions.
- ``cell_type`` – The cell type for which the differential abundance is computed.  
- ``tissue_general`` – The general tissue category associated with the dataset.  
- ``disease`` – The disease condition involved in the comparison.  
- ``baseline`` – The reference condition used for comparison (e.g., "normal").  
- ``ncell_disease`` – The number of cells sampled in the disease condition.  
- ``ncell_baseline`` – The number of cells sampled in the baseline (normal) condition.  
- ``frac_baseline`` – The proportion of the cell type in the baseline (normal) condition.  
- ``frac_disease`` – The proportion of the cell type in the disease condition.  
- ``delta_frac`` – The difference in cell type proportion between disease and baseline. 


Differential gene expression
++++++++++++++++++++++++++++++
**Endpoint**: ``/differential_gene_expression``

**Method**: ``POST``  

**Description**:

This endpoint retrieves differentially expressed genes between a baseline condition and a specified state (e.g., disease vs. normal). By default, it identifies the **top 10 up and down-regulated genes** in each cell type across all datasets that match the filter criteria.

**Parameters**:

- ``differential_axis`` *(default: disease)* – The axis along which differential expression is calculated, i.e. disease = disease condition vs normal condition.
- ``feature`` *(optional)* – The gene to query.
- ``top_n`` *(optional,default: 10)* – Number of top differentially expressed genes to return.
- ``method`` *(default: delta_fraction)* – Method to calculate differential expression (``delta_fraction`` or ``ratio_average``).
- ``disease`` *(optional)* – Filter by disease.
- ``cell_type`` *(optional)* – Filter by cell type.
- ``tissue`` *(optional)* – Filter by tissue type.
- ``sex`` *(optional)* – Filter by sex.
- ``development_stage`` *(optional)* – Filter by developmental stage.
  
.. note::

   - If ``feature`` is provided, ``top_n`` is ignored.

**Returns: A list of dictionaries, each containing the following key-value pairs:**

- ``tissue`` – The tissue where the cell was extracted from.
- ``cell_type`` – The specific cell type.
- ``regulation`` – Indicates whether the gene is **up** or **down** regulated.
- ``gene`` – The queried gene or a top-ranked differentially expressed gene.
- ``unit`` – Measurement unit (e.g., "cptt").
- ``baseline_expr`` – The gene's expression level in the baseline condition.
- ``state_expr`` – The gene's expression level in the specified condition.
- ``baseline_fraction`` – Fraction of cells expressing the gene in the baseline.
- ``state_fraction`` – Fraction of cells expressing the gene in the condition.
- ``metric`` – The computed differential expression value.
- ``dataset_id`` – The dataset from which this differential expression result was computed.
- ``differential_axis`` – The comparison category (e.g., disease, sex).
- ``state`` – The queried condition (e.g., "COVID-19").
- ``baseline`` – The reference condition (e.g., "normal").


Highest measurement
++++++++++++++++++++++++++++++
**Endpoint**: ``/highest_measurement``

**Method**: ``POST``  

**Description**:

Retrieve the top N cell types and tissue combination with the highest expression of a given feature (gene) across multiple datasets. This helps identify the most highly expressing cell types for a gene of interest in different diseases and tissues.

**Parameters**:

- ``feature`` *(required)* – The gene to query.  
- ``number`` *(optional, default: 10)* – Number of highest expressing cell types to return.  

**Returns**:  
A list of top-expressing cell types for the specified gene, ordered by expression level.

Each object contains:

- ``dataset_id`` – The dataset identifier.  
- ``cell_type`` – The cell type with high expression of the specified gene.  
- ``tissue_general`` – The broad anatomical location of the extracted cells.  
- ``disease`` – The associated disease condition (e.g., `"COVID-19"` or `"normal"`).  
- ``cell_count`` – The number of cells of this type in the dataset.  
- ``expression`` – The average expression level of the queried gene in this cell type.  

.. note::

   - The results rank the highest expressors of the queried gene based on cell type and tissue.
   - If ``number`` is greater than the available results, all possible results are returned.
  
Average
++++++++++++++++++++
**Endpoint**: ``/average``

**Method**: ``POST``  

**Description**:

Retrieve the average expression levels of one or more selected features (e.g., genes) across cell types, tissues, and diseases. This endpoint aggregates gene expression values from multiple datasets to provide an overview of average expression.

**Parameters**:

- ``features`` *(required)* – A comma-separated list of features (genes) to query.
- ``disease`` *(optional)* – Filter by disease.
- ``cell_type`` *(optional)* – Filter by cell type.
- ``tissue`` *(optional)* – Filter by tissue.
- ``sex`` *(optional)* – Filter by sex.
- ``development_stage`` *(optional)* – Filter by developmental stage.
- ``unique_ids`` *(optional)* – The unique_ids user picked from metadata result.
- ``include_normal`` *(optional, default: False)*

**Returns: A list of dictionaries, each containing:**

- ``cell_count`` – The number of cells in the given category.  
- ``cell_type`` – The cell type associated with the measurement.  
- ``tissue_general`` – The general tissue where the cell type is found.  
- ``disease`` – The disease condition associated with the measurement.  
- ``dataset_id`` – The dataset from which the measurement originates.  

.. note::

   - If ``include_normal=True``, results pair each disease entry with its corresponding normal condition, appearing consecutively.
   - Only applicable when a disease filter is provided.   
   - If no disease is specified, results naturally include both disease and normal conditions, making ``include_normal`` redundant.


Fraction detected
++++++++++++++++++
**Endpoint**: ``/fraction_detected``

**Method**: ``POST``  

**Description**:

Retrieve the fraction of cells in which a given gene is detected across different cell types, tissues, and diseases. This provides an estimation of how commonly a gene is expressed in a given cell population.

**Parameters**:

- ``features`` *(required)* – A comma-separated list of features (genes) to query.
- ``disease`` *(optional)* – Filter by disease.
- ``cell_type`` *(optional)* – Filter by cell type.
- ``tissue`` *(optional)* – Filter by tissue.
- ``sex`` *(optional)* – Filter by sex.
- ``development_stage`` *(optional)* – Filter by developmental stage.
- ``unique_ids`` *(optional)* – The unique_ids user picked from metadata result.
- ``include_normal`` *(optional, default: False)*

**Returns: A list of dictionaries, each containing:**  

- ``cell_count`` – The number of cells in the given category.  
- ``cell_type`` – The cell type associated with the measurement.  
- ``tissue_general`` – The general tissue where the cell type is found.  
- ``disease`` – The disease condition associated with the measurement.  
- ``dataset_id`` – The dataset from which the measurement originates.  
- One or more feature-specific values representing the **fraction detected** for the queried genes.

**Example Request**:
``https://api-disease.atlasapprox.org/v1/fraction_detected?disease=Covid&features=COL1A1,CXCL1,IL6``


Dot plot
+++++++++++++
**Endpoint**: ``/dotplot``

**Method**: ``POST``  

**Description**:

Retrieve both the **average expression** and **fraction detected** for a list of genes across different cell types, tissues, and diseases. This endpoint is used for visualizing gene expression in a **dot plot format**, where dot size represents fraction detected and color represents average expression.

**Parameters**:

- ``features`` *(required)* – A comma-separated list of features (genes) to query.
- ``disease`` *(optional)* – Filter by disease.
- ``cell_type`` *(optional)* – Filter by cell type.
- ``tissue`` *(optional)* – Filter by tissue.
- ``sex`` *(optional)* – Filter by sex.
- ``development_stage`` *(optional)* – Filter by developmental stage.
- ``unique_ids`` *(optional)* – The unique_ids user picked from metadata result.
- ``include_normal`` *(optional, default: False)* –  

**Returns: A list of dictionaries, each containing:**

- ``cell_count`` – The number of cells in the given category.  
- ``cell_type`` – The cell type associated with the measurement.  
- ``tissue_general`` – The general tissue where the cell type is found.  
- ``disease`` – The disease condition associated with the measurement.  
- ``dataset_id`` – The dataset from which the measurement originates.  
- Feature-specific values:
  - ``feature`` – The gene queried.
  - ``fraction_feature`` – The fraction of cells expressing the gene.
  - ``average_feature`` – The average expression of the gene.