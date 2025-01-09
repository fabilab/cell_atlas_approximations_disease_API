import requests

from atlasapprox_disease import BadRequestError

def _fetch_metadata(api, disease: str, cell_type: str, tissue: str, sex: str, development_stage_general: str):
    """
        Fetch metadata from the API.
    """
    response = requests.get(
        api.baseurl + "metadata",
        params={
            "disease": disease,
            "cell_type": cell_type,
            "tissue": tissue,
            "sex": sex,
            "development_stage_general": development_stage_general,
        }
    )
    
    if response.ok:
        return response.json()
    else:
        raise BadRequestError(response.json()["message"])