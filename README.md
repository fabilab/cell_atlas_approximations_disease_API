<img src="https://raw.githubusercontent.com/fabilab/cell_atlas_approximations/main/figures/figure_disease_API.png" width="150" height="150">

# Disease cell atlas approximations - API
This repository provides an API to fetch and analyze cell atlas data on various human diseases. For the time being, the API uses data from [Cellxgene Census](https://chanzuckerberg.github.io/cellxgene-census/). Single cell omic data is compressed into h5 files using [scquill](https://github.com/fabilab/scquill).

---
**NOTE**

Compared to Cellxgene itself, this API is designed to be much faster and closer to typical biomedical questions.

---

## Features
Currently implemented:
- Analyze differential cell type abundance across disease states, sexes, and developmental stage.
- Analyze differential gene expression for all cell types or specific cell types in disease contexts.
- Fetch metadata related to diseases, including datasets, organs, and conditions.
- Retrieve average gene expression for specific genes across disease states.
- Identify top differentially expressed genes in disease-related contexts (e.g., kidney diseases).


## Documentation
Tutorial and reference documentation is available at [https://cell-atlas-approximations-disease-api.readthedocs.io/en/latest/](https://cell-atlas-approximations-disease-api.readthedocs.io/en/latest/).


## Usage

<details> 
    <summary> REST </summary>

### REST
The REST interface is language-agnostic and can be queried using any HTTP request handler, e.g. in JavaScript:

```javascript
(async () => {
  let response = await fetch("https://api-disease.atlasapprox.org/v1/metadata?disease=covid&cell_type=alveolar type 2");
  if (response.ok) {
    let data = await response.json();
    console.log(data);
  }  
})();
```

</details>

<details>
  <summary>Python</summary>

### Python
The Python interface uses a central `API` class. Its methods implement the REST endpoints:

```python
import atlasapprox_disease

api = atlasapprox_disease.API()
print(api.metadata())
print(api.average(disease="covid", features="IGHG1,CXCL13,S100A8"))
```
</details>

<details>
  <summary>JavaScript</summary>

### JavaScript/nodejs
An object containing one function for each API endpoint is exported by the `atlasapprox` npm package:

```javascript
let atlasapprox = require('atlasapprox-disease');
(async () => {
  let data = await api.metadata(disease:"covid");
  console.log(data);
  }  
)();

```
</details>