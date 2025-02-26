<img src="https://raw.githubusercontent.com/fabilab/cell_atlas_approximations/main/figures/figure_disease_API.png" width="150" height="150">

# Cell Atlas Approximations - Disease API
This repository provides an API to fetch and analyze cell atlas data on various human diseases. For the time being, the API uses data from [Cellxgene Census](https://chanzuckerberg.github.io/cellxgene-census/). Single cell omic data is compressed into h5 files using [scquill](https://github.com/fabilab/scquill).

---
**NOTE**

Compared to Cellxgene itself, this API is designed to be much faster and closer to typical biomedical questions.

---

## Features
Currently implemented:
- Analyze differential cell type abundance.
- Analyze differential gene expression for all cell types or specific cell types.
- Fetch metadata related to diseases.

## Python package

At the moment, the preferred way to access this API is through our Python package [atlasapprox_disease](https://test.pypi.org/project/atlasapprox_disease/):

```bash
pip install atlasapprox_disease
```

For detailed examples and tutorials, please check the `tutorials` folder.
