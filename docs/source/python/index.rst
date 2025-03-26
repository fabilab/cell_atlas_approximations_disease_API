AtlasApprox Disease - Python API
================================

This Python package provides a convenient interface to access the **atlasapprox-disease** REST API, enabling efficient querying of disease-related cell atlas approximations.

Installation
------------
You can use `pip` to install the `atlasapprox-disease` package:

.. code-block:: bash

    pip install atlasapprox-disease

Getting Started
---------------
To use the API, import the `atlasapprox_disease` package and create an API instance.

.. code-block:: python

    import atlasapprox_disease as aad
    api = aad.API()

Here’s an example of querying metadata for datasets related to COVID-19:

.. code-block:: python

    metadata = api.metadata(disease="COVID")
    print(metadata.head())

Reference API
-------------
.. autoclass:: atlasapprox_disease.API
   :members:
   :undoc-members:
   :show-inheritance: