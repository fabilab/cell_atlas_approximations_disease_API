import numpy as np
import pandas as pd
import scquill

from models.metadata import get_metadata
from models.paths import get_dataset_path

from models.exceptions import (
    FeatureNotFoundError,
)

def get_highest_measurement(feature, number=10, groupby=None, **filters):
    """Compute the highest expressors of a given feature (gene) across all diseases and datasets.

    Parameters:
        feature (str): The feature of interest.
        number (int): The number of highest expressors to return.

    Returns:
        list: A list of dictionaries containing the highest expressors of the feature.
    """
    
    if type(number) == str:
        try:
            number = int(number)
        except:
            print(f"Requested number '{number}' is not valid, using default 10")
            number = 10
            
    if groupby is None:
        groupby = []
    if "cell_type" not in groupby:
        groupby = ["cell_type"] + groupby
        for i, name in enumerate(groupby):
            if name == "tissue":
                groupby[i] = "tissue_general"

    meta = get_metadata(**filters)

    result = []
    for dataset_id, obs in meta.groupby("dataset_id"):
        # All the entries in obs are guaranteed to have nonzero cells for both normal and disease
        approx = scquill.Approximation.read_h5(get_dataset_path(dataset_id))
        # NOTE:: we should binarize consistently
        adata = approx.to_anndata(
            groupby=groupby,
        )
        
        # Check if the provided feature name is valid.
        if feature not in adata.var_names:
            raise FeatureNotFoundError(
                msg=f"Feature '{feature}' not found in dataset '{dataset_id}'.",
                feature=feature
            )

        obs_names = obs[groupby].agg("\t".join, axis=1).values
        obs_names = pd.Index(obs_names).drop_duplicates()
        adata.obs["dataset_id"] = dataset_id
        adata.obs["expression"] = np.asarray(adata[:, feature].X).ravel()
        if number > adata.n_obs:
            res = adata.obs
        else:
            res = adata.obs.nlargest(number, "expression")
        result.append(res)

    result = pd.concat(result)
    result = result.nlargest(number, "expression")

    return result
