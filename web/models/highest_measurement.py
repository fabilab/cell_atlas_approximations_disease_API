import numpy as np
import pandas as pd
import scquill

from models.metadata import get_metadata
from models.paths import get_dataset_path


def get_highest_measurement(feature, number=10, **filters):
    """Compute the highest expressors of a given feature (gene) across all diseases and datasets.

    Parameters:
        feature (str): The feature of interest.
        number (int): The number of highest expressors to return.

    Returns:
        list: A list of dictionaries containing the highest expressors of the feature.
    """

    meta = get_metadata(**filters)

    result = []
    for dataset_id, obs in meta.groupby("dataset_id"):
        # All the entries in obs are guaranteed to have nonzero cells for both normal and disease
        approx = scquill.Approximation.read_h5(get_dataset_path(dataset_id))
        # NOTE:: we should binarize consistently
        adata = approx.to_anndata(
            groupby=["tissue_general", "cell_type", "disease"],
        )
        obs_names = (
            obs[["tissue_general", "cell_type", "disease"]]
            .agg("\t".join, axis=1)
            .values
        )
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
