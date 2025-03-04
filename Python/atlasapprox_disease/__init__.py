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
    """Main object used to access the disease approximation API"""

    cache = {}

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
    ) -> pd.DataFrame:
        """Fetch metadata based on various filters

        Args:
            disease: Filter by disease name (e.g., "flu", optional).
            cell_type: Filter by cell type (optional).
            tissue: Filter by tissue (optional).
            sex: Filter by sex (e.g., "male" or "female", optional).
            development_stage: Filter by development stage (e.g., 'adult', optional).

        Returns:
            pd.DataFrame: A DataFrame containing the metadata.
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
    ) -> pd.DataFrame:
        """Get differential cell type abundance between conditions.

        Args:
            differential_axis: The axis to compute differential abundance on (default: "disease")
            disease: Filter by disease name (optional)
            cell_type: Filter by cell type (optional)
            tissue: Filter by tissue (optional)
            sex: Filter by sex (optional)
            development_stage: Filter by development stage (optional)

        Returns:
            pd.DataFrame: A DataFrame containing the differential cell type abundance.
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
    ) -> pd.DataFrame:
        """Get differential gene expression between two conditions.

        Args:
            differential_axis: The axis to compute differential abundance on (default: "disease")
            disease: Filter by disease name (optional)
            cell_type: Filter by cell type (optional)
            tissue: Filter by tissue (optional)
            sex: Filter by sex (optional)
            development_stage: Filter by development stage (optional)
            top_n: Top N differentially UP regulated genes +  Top N differentially DOWN regulated genes  (default: 10)
            feature: Specific feature to query. (optional)
            method: Method of calculation ('delta_fraction' | 'ratio_average').
        
        Returns:
            pd.DataFrame: A DataFrame containing differential gene expression.
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
    )-> pd.DataFrame:
        """
        Get the highest measurement of a specific feature.

        Args:
            feature (str): The feature (gene) to query.
            number (int): The number of highest expressors to return.

        Returns:
            pd.DataFrame: A DataFrame containing the highest measurements.
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
    ) -> pd.DataFrame:
        """
        Get the average expression of given genes.

        Args:
            features (str]): The features (genes) to query.
            disease (str): Filter by disease name.
            cell_type (str): Filter by cell type.
            tissue (str): Filter by tissue.
            sex (str): Filter by sex.
            development_stage (str): Filter by development stage.
            unique_ids (str): Filter by unique_ids from metadata.
            include_normal (bool): Include normal condition when querying a disease.

        Returns:
            pd.DataFrame: A DataFrame containing the average expression.
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
    ) -> pd.DataFrame:
        """
        Get the fraction of a given gene detected in datasets.

        Args:
            features (str]): The features (genes) to query.
            disease (str): Filter by disease name.
            cell_type (str): Filter by cell type.
            tissue (str): Filter by tissue.
            sex (str): Filter by sex.
            development_stage (str): Filter by development stage.
            unique_ids (str]): Filter by unique_ids from metadata.
            include_normal (bool): Include normal condition when querying a disease.

        Returns:
            pd.DataFrame: A DataFrame containing the fraction detected.
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
    ) -> pd.DataFrame:
        """
        Prepare data for a dotplot including average expression and fraction detected.

        Args:
            features (str): The features (genes) to query.
            disease (str): Filter by disease name.
            cell_type (str): Filter by cell type.
            tissue (str): Filter by tissue.
            sex (str): Filter by sex.
            development_stage (str): Filter by development stage.
            unique_ids (str): Filter by specific dataset IDs.
            include_normal (bool): Include normal condition when querying a disease.

        Returns:
            pd.DataFrame: A DataFrame containing dotplot data.
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