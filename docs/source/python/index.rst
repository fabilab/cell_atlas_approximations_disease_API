Python
================================

The Python interface can be used to access the disease atlas approximation API from Python. It enables efficient querying of disease-related cell atlas approximations.

Requirements
------------
You need the following Python packages:

  - ``requests``
  - ``pandas``

Installation
------------
You can use `pip` to install the `atlasapprox_disease` package:

.. code-block:: bash

    pip install atlasapprox-disease

Getting Started
---------------
To use the API, import the `atlasapprox_disease` package:

.. code-block:: python

    import atlasapprox_disease as aad

and initialise the ``API`` object:

.. code-block:: python
    
    api = aad.API()

Here’s an example of querying metadata for datasets related to COVID-19 in lung tissue and then using the `unique_ids` to query average gene expression:

.. code-block:: python

    # Step 1: Query metadata to get unique_ids
    metadata = api.metadata(disease="covid", tissue="lung")
    print(metadata.head())

    # Step 2: Use a unique_id to query average expression of specific genes
    unique_id = metadata["unique_id"].iloc[0]  # Select the first unique_id
    avg_expr = api.average(features="IGHG1,CXCL13,S100A8", unique_ids=unique_id)
    print(avg_expr)

.. note::

   When using `unique_ids` in methods like `average`, `fraction_detected`, or `dotplot`, only specify the `features` parameter alongside it. Do not include other metadata filters (`disease`, `cell_type`, `tissue`, `sex`, `development_stage`), as `unique_ids` already encapsulate these conditions. Combining them will raise a `ParamsConflictError`.

Examples
--------
Explore practical examples of using the Python API to analyze disease-related single-cell data:

.. toctree::
   :maxdepth: 1

   Python Examples <../auto_examples/index>

Reference API
-------------
.. autoclass:: atlasapprox_disease.API
   :members:
   :undoc-members:
   :show-inheritance: