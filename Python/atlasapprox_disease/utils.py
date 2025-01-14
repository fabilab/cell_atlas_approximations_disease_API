import pandas as pd
import requests
from typing import Union, List

from atlasapprox_disease import BadRequestError


def _fetch_metadata(
    api, disease: str, cell_type: str, tissue: str, sex: str, development_stage: str
):
    """
    Fetch metadata from the API.
    
    Return:
        A pandas.DataFrame with metadata that satisfy the filters
    """
    response = requests.get(
        api.baseurl + "metadata",
        params={
            "disease": disease,
            "cell_type": cell_type,
            "tissue": tissue,
            "sex": sex,
            "development_stage": development_stage,
        },
    )

    if response.ok:
        return response.json()
    else:
        raise BadRequestError(response.json()["message"])


def _fetch_differential_cell_type_abundance(
    api,
    differential_axis: str,
    disease: str,
    cell_type: str,
    tissue: str,
    sex: str,
    development_stage: str,
    unique_ids: Union[str, List[str]],
):
    """
    Fetch differential cell type abundance data from the API.
    
    Returns:
        A pandas.DataFrame with the differential cell type abundance results.
    """
    response = requests.post(
        api.baseurl + "differential_cell_type_abundance",
        params={
            "differential_axis": differential_axis,
            "disease": disease,
            "cell_type": cell_type,
            "tissue": tissue,
            "sex": sex,
            "development_stage": development_stage,
            "unique_ids": unique_ids,
        },
    )
    if response.ok:
        resjson = response.json()
        df = pd.DataFrame(resjson)
        return df
    else:
        raise BadRequestError(response.json()["message"])


def _fetch_differential_gene_expression(
    api,
    differential_axis: str,
    disease: str,
    cell_type: str,
    tissue: str,
    sex: str,
    development_stage: str,
    top_n: int,
    feature: str,
    method: str,
    unique_ids: Union[str, List[str]],
):
    """
    Fetch top N differential expressed genes or expression of queried genes from the API.
    
    Returns:
        A pandas.DataFrame with differential gene expression results.

    """
    response = requests.post(
        api.baseurl + "differential_gene_expression",
        params={
            "differential_axis": differential_axis,
            "disease": disease,
            "cell_type": cell_type,
            "tissue": tissue,
            "sex": sex,
            "development_stage": development_stage,
            "top_n": top_n,
            "feature": feature,
            "method": method,
            "unique_ids": unique_ids,
        },
    )
    if response.ok:
        resjson = response.json()
        df = pd.DataFrame(resjson)
        return df
    else:
        raise BadRequestError(response.json()["message"])
