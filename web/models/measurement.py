import numpy as np
import pandas as pd
import scquill

from models.metadata import get_metadata
from models.paths import get_dataset_path
from models.exceptions import SomeFeaturesNotFoundError


def get_measurement(features, kind, groupby=None, **filters):
    """Compute the highest expressors of a given feature (gene) across all diseases and datasets.

    Parameters:
        features (list): List of genes/features to query.
        kind (string): The type of measurement (e.g., "average", "fraction").
        groupby (optional): Variables to group by when aggregating.
        filters (dict): Metadata filters, including `unique_ids` if provided.

    Returns:
        pd.DataFrame
    """
    if groupby is None:
        groupby = []
        
    # Track whether unique_ids were provided
    use_unique_ids = "unique_ids" in filters

    if use_unique_ids:
        # remove empty space
        unique_ids = [uid.strip() for uid in filters.pop("unique_ids") if isinstance(uid, str)]

        # Retrieve metadata for the given `unique_ids`
        metadata = get_metadata()
        meta = metadata[metadata["unique_id"].isin(unique_ids)]

        if meta.empty:
            raise ValueError(f"No metadata found for the provided unique IDs: {unique_ids}")

        # Create a mapping of metadata per unique_id
        metadata_map = meta.set_index("unique_id").to_dict(orient="index")
        # Extract dataset-specific filters from metadata
        metadata_columns = [col for col in meta.columns if col not in {"unique_id", "dataset_id", "cell_count"}]
        filters = {col: meta[col].unique().tolist() for col in metadata_columns}
    
    else:
        meta = get_metadata(**filters)
        metadata_map = {} 
        metadata_columns = list(filters.keys())
  
    if "cell_type" not in groupby:
        groupby = ["cell_type"] + groupby
        for i, name in enumerate(groupby):
            if name == "tissue":
                groupby[i] = "tissue_general"

    # For each filter, we need to make sure that the data is not coarse grained before we apply it
    filter_cols = list(filters.keys())
    filter_cols_additional = [col for col in filter_cols if col not in groupby]
    groupby_with_filters = groupby + filter_cols_additional
    
    # Determine whether to group by `unique_id` or `dataset_id`
    grouping_column = "unique_id" if use_unique_ids else "dataset_id"
    
    result = []
    for group_key, obs in meta.groupby(grouping_column):
        # Retrieve the correct dataset_id (important when grouping by `unique_id`)
        dataset_id = obs["dataset_id"].values[0]

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

        # Ensure `dataset_id` is included in the response
        adata.obs["dataset_id"] = dataset_id

        # Filter and coarse grain, in that order
        if len(filter_cols) > 0 and filter_cols != ["dataset_id"]:
            obs_filter_unique = obs.set_index(filter_cols).index.drop_duplicates()
            idx_obs = adata.obs_names[
                adata.obs.set_index(filter_cols).index.isin(obs_filter_unique)
            ]
            adata = adata[idx_obs]
            adata = scquill.utils.coarse_grain_anndata(adata, groupby=groupby)
            
        res = adata.obs.copy()

        # Assign metadata per `unique_id`
        if use_unique_ids:
            if group_key in metadata_map:
                for field, value in metadata_map[group_key].items():
                    res.loc[:, field] = value


        # Assign metadata per dataset for non-unique_id queries
        else:
            matching_meta = meta[meta["dataset_id"] == dataset_id]
            for field in metadata_columns + ["dataset_id"]:  # Ensure dataset_id is included
                if field in matching_meta.columns and field not in res.columns:
                    res[field] = matching_meta[field].values[0]

        # Add gene expression values
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
