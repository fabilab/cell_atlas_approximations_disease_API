import pandas as pd
import requests

def _fetch_metadata(api, **kwargs):
    """
    Fetch metadata from the API.

    Return:
        A pandas.DataFrame with metadata that satisfy the filters
    
    Raises:
        BadRequestError: If the API request fails.
    """
    
    response = requests.get(
        api.baseurl + "metadata",
        params=kwargs
    )

    if response.ok:
        return pd.DataFrame(response.json())
    else:
        return response.json()


def _fetch_differential_cell_type_abundance(api, **kwargs):
    """
    Fetch differential cell type abundance data from the API.

    Returns:
        A pandas.DataFrame with the differential cell type abundance results.
    """
    response = requests.post(
        api.baseurl + "differential_cell_type_abundance",
        params=kwargs
    )
    if response.ok:
        return pd.DataFrame(response.json())
    else:
        return response.json()


def _fetch_differential_gene_expression(api, **kwargs):
    """
    Fetch top N differential expressed genes or expression of queried genes from the API.

    Returns:
        A pandas.DataFrame with differential gene expression results.

    """
    response = requests.post(
        api.baseurl + "differential_gene_expression",
        params=kwargs
    )
    if response.ok:
        return pd.DataFrame(response.json())
    else:
        return response.json()



def _fetch_highest_measurement(api, **kwargs):
    """
    Fetch the highest measurement of a specific feature from the API.

    Return:
        A DataFrame containing the highest measurements.
    """
    response = requests.post(
        api.baseurl + "highest_measurement",
        params=kwargs
    )
    if response.ok:
        return pd.DataFrame(response.json())
    else:
        return response.json()

def _fetch_average(api, **kwargs):
    """
    Fetch the average expression of specific features from the API.
    
    Returns:
        A DataFrame containing the average expression values.
    """
    response = requests.post(
        api.baseurl + "average",
        params=kwargs
    )
    if response.ok:
        return pd.DataFrame(response.json())
    else:
        return response.json()
    
def _fetch_fraction_detected(api, **kwargs):
    """Fetch the fraction of specific features detected in datasets.

    Args:
        api: The API object making the request.
        **kwargs: Optional filter arguments.

    Returns:
        A DataFrame containing the fraction detected.
    """
    response = requests.post(
        api.baseurl + "fraction_detected",
        params=kwargs
    )
    if response.ok:
        return pd.DataFrame(response.json())
    else:
        return response.json()

def _fetch_dotplot(api, **kwargs):
    """Fetch dot plot data for the specified features from the API.

    Args:
        api: The API object making the request.
        **kwargs: Optional filter arguments.

    Returns:
        A DataFrame containing dotplot data.
    """
        
    response = requests.post(
        api.baseurl + "dotplot",
        params=kwargs
    )
    if response.ok:
        return pd.DataFrame(response.json())
    else:
        return response.json()