import numpy as np
import pandas as pd
import scquill

from models.metadata import get_metadata
from models.paths import get_dataset_path
from models.exceptions import SomeFeaturesNotFoundError


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

    # For each filter, we need to make sure that the data is not coarse grained before we apply it
    filter_cols = list(filters.keys())
    filter_cols_additional = [col for col in filter_cols if col not in groupby]
    groupby_with_filters = groupby + filter_cols_additional

    result = []
    for dataset_id, obs in meta.groupby("dataset_id"):
        # All the entries in obs are guaranteed to have nonzero cells for both normal and disease
        approx = scquill.Approximation.read_h5(get_dataset_path(dataset_id))

        # NOTE:: we should binarize consistently
        adata = approx.to_anndata(
            groupby=groupby_with_filters,
        )

        # Validate feature names.
        # This validation slows down the request by approximately 2 seconds. Optimization needed in the future.
        invalid_features = [feature for feature in features if feature not in adata.var_names]
        if invalid_features:
            raise SomeFeaturesNotFoundError(
                msg=f"Some features could not be found in dataset '{dataset_id}': {', '.join(invalid_features)}.",
                features=invalid_features
            )
        
        adata = approx.to_anndata(
            groupby=groupby_with_filters,
            features=features
        )
            
        adata.obs["dataset_id"] = dataset_id

        # Filter and coarse grain, in that order
        # dataset_id is intrinsically filtered by the for loop
        if len(filter_cols) > 0 and filter_cols != ["dataset_id"]:
            obs_filter_unique = obs.set_index(filter_cols).index.drop_duplicates()
            idx_obs = adata.obs_names[
                adata.obs.set_index(filter_cols).index.isin(obs_filter_unique)
            ]
            adata = adata[idx_obs]
            adata = scquill.utils.coarse_grain_anndata(adata, groupby=groupby)

        # Now we can index the adata properly
        obs_names = obs[groupby].agg("\t".join, axis=1).values
        obs_names = pd.Index(obs_names).drop_duplicates()

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
