"""
Cell atlas approximations(disease) - Python API Interface
"""

import os
import pandas as pd
from typing import Union, List

from atlasapprox_disease.exceptions import BadRequestError
from atlasapprox_disease.utils import (
    _fetch_metadata,
    _fetch_differential_cell_type_abundance,
    _fetch_differential_gene_expression,
    _fetch_highest_measurement,
    _fetch_average,
    _fetch_fraction_detected,
    _fetch_dotplot,
)


__version__ = "0.1.4"

__all__ = (
    "api_version",
    "API",
    "BadRequestError",
    __version__,
)
print("Loading atlasapprox_disease from:", __file__)
api_version = "v1"

baseurl = os.getenv(
    "ATLASAPPROX_DISEASE_BASEURL",
    "https://api-disease.atlasapprox.org",
)
baseurl = baseurl.rstrip("/") + "/"
baseurl += f"{api_version}/"

show_credit = os.getenv("ATLASAPPROX_DISEASE_HIDECREDITS") is None
credit = """Data sources for the disease approximations:
    CellxGene Census (https://chanzuckerberg.github.io/cellxgene-census/)

To hide this message, set the environment variable ATLASAPPROX_DISEASE_HIDECREDITS to any
nonzero value, e.g.:

import os
os.environ[ATLASAPPROX_DISEASE_HIDECREDITS"] = "yes"
import atlasapprox_disease

To propose a new disease be added to the list of approximations, please contact
Fabio Zanini (fabio DOT zanini AT unsw DOT edu DOT au).
"""

if show_credit:
    print(credit)
    show_credit = False


class API:
    """Main object used to access the atlasapprox-disease REST API."""

    def __init__(self, url=None):
        """Create an instance of the atlasapprox_disease API."""
        self.baseurl = url if url is not None else baseurl

    def metadata(
        self,
        disease: str = None,
        cell_type: str = None,
        tissue: str = None,
        sex: str = None,
        development_stage: str = None,
    ):
        """Retrieves metadata records from the atlasapprox-disease API. Each record
        represents a unique combination of dataset, cell type, tissue, disease condition, sex, and
        developmental stage that meets the query criteria.
        
        Args:
            disease (str, optional): Filter by disease name (e.g., "covid").
            cell_type (str, optional): Filter by cell type (e.g., "fibroblast").
            tissue (str, optional): Filter by tissue (e.g., "lung"). 
            sex (str, optional): Filter by sex (e.g., "male", "female").
            development_stage (str, optional): Filter by developmental stage (e.g., "adult").

        Returns:
            A DataFrame containing metadata records for datasets matching the filters.
        """
        return _fetch_metadata(self, disease=disease, cell_type=cell_type, tissue=tissue, sex=sex, development_stage=development_stage)

    def differential_cell_type_abundance(
        self,
        differential_axis: str = None,
        disease: str = None,
        cell_type: str = None,
        tissue: str = None,
        sex: str = None,
        development_stage: str = None,
    ):
        """Get differential cell type abundance between a specified condition (e.g., disease)
        and a baseline (e.g., normal) across datasets.

        Args:
            differential_axis (str, optional): Axis for comparison (default: "disease"). Options: "disease" (disease vs. normal), "sex" (male vs. female).
            disease (str, optional): Filter by disease name (e.g., "flu").
            cell_type (str, optional): Filter by cell type (e.g., "macrophage").
            tissue (str, optional): Filter by tissue (e.g., "lung").
            sex (str, optional): Filter by sex (e.g., "male", "female").
            development_stage (str, optional): Filter by developmental stage (e.g., "adult").

        Returns:
            A DataFrame containing differential cell type abundance between conditions.
        """
        return _fetch_differential_cell_type_abundance(
            self,
            differential_axis=differential_axis,
            disease=disease,
            cell_type=cell_type,
            tissue=tissue,
            sex=sex,
            development_stage=development_stage,
        )

    def differential_gene_expression(
        self,
        differential_axis: str = None,
        disease: str = None,
        cell_type: str = None,
        tissue: str = None,
        sex: str = None,
        development_stage: str = None,
        top_n: int = None,
        feature: str = None,
        method: str = None,
    ):
        """Get differential gene expression between conditions.

        Args:
            differential_axis (str, optional): Axis for comparison (default: "disease").
                Options: "disease" (disease vs. normal), "sex" (male vs. female), "age" (e.g., adult vs. other).
            disease (str, optional): Filter by disease name (e.g., "COVID").
            cell_type (str, optional): Filter by cell type (e.g., "T cell").
            tissue (str, optional): Filter by tissue (e.g., "lung").
            sex (str, optional): Filter by sex (e.g., "male", "female").
            development_stage (str, optional): Filter by developmental stage (e.g., "adult").
            
            top_n (int, optional): Number of top up- and down-regulated genes to return (default: 10).
                Ignored if `feature` is provided.
            feature (str, optional): Specific gene to query (e.g., "IL6"). If provided, `top_n` is ignored.
            method (str, optional): Method to compute differential expression (default: "delta_fraction").
                Options: "delta_fraction", "ratio_average".
        
        Returns:
            A DataFrame containing differential gene expression results between conditions.
        """
        return _fetch_differential_gene_expression(
            self,
            differential_axis=differential_axis,
            disease=disease,
            cell_type=cell_type,
            tissue=tissue,
            sex=sex,
            development_stage=development_stage,
            top_n=top_n,
            feature=feature,
            method=method,
        )
        
    def highest_measurement(
        self,
        feature: str =  None,
        number : int = None,
     ) :
        """
        Retrieves the top N cell types and tissue combinations with the highest expression
        of a specified gene across datasets. It helps identify which cell types most highly express a
        gene of interest in different diseases and tissues.

        Args:
            feature (str): The gene to query (e.g., "IL6").
            number (int, optional): Number of top-expressing cell types to return (default: 10).

        Returns:
            A DataFrame containing the top cell types with the highest expression of the specified feature.
        """
        return _fetch_highest_measurement(self, feature=feature, number=number)

    def average(
        self,
        features: str = None,
        disease: str = None,
        cell_type: str = None,
        tissue: str = None,
        sex: str = None,
        development_stage: str = None,
        unique_ids: str = None,
        include_normal: bool = None
    ):
        """
        Get the average expression levels of one or more genes across cell types,
        tissues, sex, development stage and diseases.

        Args:
            features (str): A comma-separated list of genes to query (e.g., "IL6,AGT").
            disease (str, optional): Filter by disease name (e.g., "diabete").
            cell_type (str, optional): Filter by cell type (e.g., "T cell").
            tissue (str, optional): Filter by tissue (e.g., "lung").
            sex (str, optional): Filter by sex (e.g., "male", "female").
            development_stage (str, optional): Filter by developmental stage (e.g., "adult").
            
            unique_ids (str, optional): A comma-separated list of unique IDs from metadata results
                to filter specific dataset entries.
            include_normal (bool, optional): If True, includes the corresponding normal condition
                alongside the queried disease (default: False). Only applicable when a disease filter
                is provided. If no disease is specified, results include both disease and normal conditions.

        Returns:
            A DataFrame containing average expression data for the specified features across conditions.
        """
        return _fetch_average(
            self,
            features=features,
            disease=disease,
            cell_type=cell_type,
            tissue=tissue,
            sex=sex,
            development_stage=development_stage,
            unique_ids=unique_ids,
            include_normal=include_normal
        )
    
    def fraction_detected(
        self,
        features: str,
        disease: str = None,
        cell_type: str = None,
        tissue: str = None,
        sex: str = None,
        development_stage: str = None,
        unique_ids: str = None,
        include_normal: bool = None
    ):
        """
        Get the fraction of cells expressing the specified features across conditions.

        Args:
            features (str): A comma-separated list of genes to query.
            disease (str, optional): Filter by disease name.
            cell_type (str, optional): Filter by cell type.
            tissue (str, optional): Filter by tissue.
            sex (str, optional): Filter by sex.
            development_stage (str, optional): Filter by developmental stage.
            
            unique_ids (str, optional): A comma-separated list of unique IDs from metadata results
                to filter specific dataset entries.
            include_normal (bool, optional): If True, includes the corresponding normal condition
                alongside the queried disease (default: False). Only applicable when a disease filter
                is provided. If no disease is specified, results include both disease and normal conditions.

        Returns:
            A DataFrame containing the fraction of cells expressing the specified features across conditions.
        """
        return _fetch_fraction_detected(
            self,
            features=features,
            disease=disease,
            cell_type=cell_type,
            tissue=tissue,
            sex=sex,
            development_stage=development_stage,
            unique_ids=unique_ids,
            include_normal=include_normal
        )
        
    def dotplot(
        self,
        features: str,
        disease: str = None,
        cell_type: str = None,
        tissue: str = None,
        sex: str = None,
        development_stage: str = None,
        unique_ids: str = None,
        include_normal: bool = None
    ):
        """
        Pepare data for a dot plot, including average expression and fraction detected. The data is suitable for visualizing in a
        dot plot format, where dot size represents fraction detected and color represents average
        expression.

        Args:
            features (str): A comma-separated list of genes to query.
            disease (str, optional): Filter by disease name.
            cell_type (str, optional): Filter by cell type.
            tissue (str, optional): Filter by tissue.
            sex (str, optional): Filter by sex.
            development_stage (str, optional): Filter by developmental stage.
            
            unique_ids (str, optional): A comma-separated list of unique IDs from metadata results to filter specific dataset entries.
            include_normal (bool, optional): If True, includes the corresponding normal condition
                alongside the queried disease (default: False). Only applicable when a disease filter
                is provided. If no disease is specified, results include both disease and normal conditions.

        Returns:
            A DataFrame containing dot plot data with average expression and fraction detected for the specified features.
        """
        return _fetch_dotplot(
            self,
            features=features,
            disease=disease,
            cell_type=cell_type,
            tissue=tissue,
            sex=sex,
            development_stage=development_stage,
            unique_ids=unique_ids,
            include_normal=include_normal
        )