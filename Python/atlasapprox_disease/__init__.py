"""
Cell atlas approximations(disease) - Python API Interface
"""

import os

from atlasapprox_disease.exceptions import BadRequestError
from atlasapprox_disease.utils import (
    _fetch_metadata
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
        development_stage_general: str = None,
    ):
        """Fetch metadata based on various filters

        Args:
            disease: Filter by disease name (e.g., "flu", optional).
            cell_type: Filter by cell type (optional).
            tissue: Filter by tissue (optional).
            sex: Filter by sex (e.g., "male" or "female", optional).
            development_stage_general: Filter by development stage (e.g., 'adult', optional).

        Returns:
            A list of metadata records matching the filters.
        """
        return _fetch_metadata(self, disease, cell_type, tissue, sex, development_stage_general)
