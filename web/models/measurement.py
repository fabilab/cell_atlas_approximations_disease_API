import numpy as np
import pandas as pd
import scquill

from models.metadata import get_metadata
from models.paths import get_dataset_path


def get_measurement(features, kind, groupby=None, **filters):
    """Compute the highest expressors of a given feature (gene) across all diseases and datasets.

    Parameters:
        features (str): The feature of interest.

    Returns:
        pd.DataFrame
    """
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
            features=features,
        )
        obs_names = obs[groupby].agg("\t".join, axis=1).values
        obs_names = pd.Index(obs_names).drop_duplicates()
        adata.obs["dataset_id"] = dataset_id

        res = adata.obs.copy()
        for feature in features:
            if kind == "fraction":
                resi = np.asarray(adata[:, feature].layers["fraction"]).ravel()
            else:
                resi = np.asarray(adata[:, feature].X).ravel()
            res[feature] = resi
        result.append(res)

    result = pd.concat(result)

    return result


def get_average(features, groupby=None, **filters):
    return get_measurement(features, kind="average", groupby=groupby, **filters)


def get_fraction_detected(features, groupby=None, **filters):
    return get_measurement(features, kind="fraction", groupby=groupby, **filters)


def get_dotplot(features, groupby=None, **filters):
    """Get both average measurement and fraction detected at once."""
    average = get_average(
        features,
        groupby=groupby,
        **filters,
    )
    fraction = get_fraction_detected(
        features,
        groupby=groupby,
        **filters,
    )
    fraction.rename(columns={fea: f"fraction_{fea}" for fea in features}, inplace=True)
    result = pd.merge(left=average, right=fraction)
    return result
