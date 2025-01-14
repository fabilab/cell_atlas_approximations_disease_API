"""
Cell atlas approximations(disease) - Python API Interface
"""

import os
from typing import Union, List

from atlasapprox_disease.exceptions import BadRequestError
from atlasapprox_disease.utils import (
    _fetch_metadata,
    _fetch_differential_cell_type_abundance,
    _fetch_differential_gene_expression,
)

__version__ = "0.0.1"

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
    ):
        """Fetch metadata based on various filters

        Args:
            disease: Filter by disease name (e.g., "flu", optional).
            cell_type: Filter by cell type (optional).
            tissue: Filter by tissue (optional).
            sex: Filter by sex (e.g., "male" or "female", optional).
            development_stage: Filter by development stage (e.g., 'adult', optional).

        Returns:
            A list of metadata records matching the filters.
        """
        return _fetch_metadata(self, disease, cell_type, tissue, sex, development_stage)

    def differential_cell_type_abundance(
        self,
        differential_axis: str = "disease",
        disease: str = None,
        cell_type: str = None,
        tissue: str = None,
        sex: str = None,
        development_stage: str = None,
        unique_ids: Union[str, List[str]] = None,
    ):
        """Get differential cell type abundance between conditions.

        Args:
            differential_axis: The axis to compute differential abundance on (default: "disease")
            disease: Filter by disease name (optional)
            cell_type: Filter by cell type (optional)
            tissue: Filter by tissue (optional)
            sex: Filter by sex (optional)
            development_stage: Filter by development stage (optional)
            unique_ids: Filter by specific dataset IDs. Can be a comma-separated string or list of strings (optional)

        Returns:
            A list of differential abundance results

        Raises:
            ValueError: If both unique_ids and other filters are specified
        """
        # Validate that unique_ids is not used with other filters
        if unique_ids is not None and any(
            [
                x is not None
                for x in [disease, cell_type, tissue, sex, development_stage]
            ]
        ):
            raise ValueError(
                "You can specify either unique_ids or metadata filters, not both"
            )

        return _fetch_differential_cell_type_abundance(
            self,
            differential_axis,
            disease,
            cell_type,
            tissue,
            sex,
            development_stage,
            unique_ids,
        )

    def differential_gene_expression(
        self,
        differential_axis: str = "disease",
        disease: str = None,
        cell_type: str = None,
        tissue: str = None,
        sex: str = None,
        development_stage: str = None,
        top_n: int = 10,
        feature: str = None,
        method: str = "delta_fraction",
        unique_ids: Union[str, List[str]] = None,
    ):
        """Get differential gene expression between disease and normal conditions.

        Args:
            differential_axis: The axis to compute differential abundance on (default: "disease")
            disease: Filter by disease name (optional)
            cell_type: Filter by cell type (optional)
            tissue: Filter by tissue (optional)
            sex: Filter by sex (optional)
            development_stage: Filter by development stage (optional)
            top_n: Top N differentially UP regulated genes +  Top N differentially DOWN regulated genes  (default: 10)
            feature: Query expression level difference for a given feature (optional)
            method: Calculation of differential expression [delta_fraction|ratio_average] (default: delta_fraction)
            unique_ids: Filter by specific dataset IDs. Can be a comma-separated string or list of strings (optional)

        Returns:
            A list of differentially expressed genes

        Raises:
            BadRequestError: If both top_n and feature are specified
        """
        return _fetch_differential_gene_expression(
            self,
            differential_axis,
            disease,
            cell_type,
            tissue,
            sex,
            development_stage,
            top_n,
            feature,
            method,
            unique_ids,
        )
