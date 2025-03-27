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
You can use `pip` to install the `atlasapprox-disease` package:

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

Here’s an example of querying metadata for datasets related to COVID-19:

.. code-block:: python

    metadata = api.metadata(disease="covid")
    print(metadata.head())

Reference API
-------------
.. autoclass:: atlasapprox_disease.API
   :members:
   :undoc-members:
   :show-inheritance: