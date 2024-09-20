# cell_atlas_approximations_disease_API

## Overview

This repository provides an API to fetch and analyze cell atlas data on various human diseases.

The API sources compressed h5 files created with the [scquill](https://github.com/fabilab/scquill) package on data from [Cellxgene Census](https://chanzuckerberg.github.io/cellxgene-census/). Compared to Cellxgene itself, this API is designed to be much faster and closer to typical biomedical questions.

## Features
Currently implemented:
- Analyze differential cell type abundance.
- Analyze differential gene expression for all cell types or specific cell types.
- Fetch metadata related to diseases.

## Python package

At the moment, the preferred way to access this API is through our Python package [atlasapprox-disease](https://test.pypi.org/project/atlasapprox-disease/):

```bash
pip install atlasapprox-diseases
```

For detailed examples and tutorials, please check the `tutorials` folder.
