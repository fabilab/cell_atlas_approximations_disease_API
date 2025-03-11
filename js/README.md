# Cell Atlas Approximations - Disease JavaScript API

**atlasapprox-disease** is a JavaScript interface designed to provide easy access to disease-related cell atlas data. Built using the cellxgene census as the initial data source and powered by the scquill approximation algorithm, this package enables researchers to explore thousands of datasets efficiently.

This project enables biologists, doctors, and data scientist to quickly find answers for questions such as:

- *Which cell types are affected by diabetes in the pancreas?*
- *How does gene expression change in macrophages during COVID-19 infection?*
- *Which cell types show the highest expression of a given gene in a disease condition?*
- *What fraction of cells express a particular gene in disease versus normal conditions?*

## Version
The latest API version is `v1`.

This package provides access to a variety of disease-related single-cell datasets, covering multiple cell types, tissues, development stage and sex.

## Installation
To install the package via npm:

```sh
npm install @fabilab/atlasapprox-disease
```

## Usage
An object containing one function for each API endpoint is exported by the `@fabilab/atlasapprox-disease` npm package:

```javascript

import atlasapprox_disease from '@fabilab/atlasapprox-disease';

(async () => {
  let data = await atlasapprox_disease.average({ features: "PRDM1,PTPRC,ACTB,GAPDH" });
  console.log(data);
})();
```

## Available API Functions
- **`metadata()`** - Retrieve unique combinations of dataset metadata, including cell types, tissues, development stage, sex and diseases.
- **`differential_gene_expression()`** - Identify differentially expressed genes between conditions.
- **`differential_cell_type_abundance()`** - Compare cell type abundance in different disease states.
- **`highest_measurement()`** - Find the cell types with the highest expression of a given gene.
- **`average()`** - Get the average expression of genes across conditions.
- **`fraction_detected()`** - Determine the fraction of cells expressing a specific gene.
- **`dotplot()`** - Retrieve dot plot data for gene expression across cell types.

## Authors
- [Fabio Zanini @ fabilab](https://fabilab.org)
- [Ying Xu @ fabilab](https://fabilab.org/pages/people.html)

