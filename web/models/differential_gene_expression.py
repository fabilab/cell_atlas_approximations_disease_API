import numpy as np
import pandas as pd
import scquill

from models.metadata import get_metadata_with_normal
from models.paths import get_dataset_path


def get_diff_expression(number=10, feature=None, method="delta_fraction", **filters):
    """
    Computes the differential expression for a given dataset.

    Parameters:
        number (int): The number of differentially expressed features to return. Either this or "feature" must be provided.
        feature (str): The feature to compute differential expression for. Either this or "number" must be provided.
        **filters: Arbitrary keyword arguments to filter the metadata.

    Returns:
        result (list): List of dictionaries containing the computed results.
    """
    if (number is not None and number > 0) and (feature is not None):
        raise ValueError("Either 'number' or 'feature' must be provided, not both.")
    if (number is None or number <= 0) and (feature is None):
        raise ValueError("Either 'number' or 'feature' must be provided.")

    meta = get_metadata_with_normal(**filters)
    if len(meta) == 0:
        return pd.DataFrame()

    # NOTE: this should be improved
    meta["disease_bin"] = meta["disease"] != "normal"

    ncells = (
        meta.groupby(["dataset_id", "tissue_general", "cell_type", "disease_bin"])[
            "cell_count"
        ]
        .sum()
        .unstack("disease_bin", fill_value=0)
    )
    # If normal is missing, we cannot compute differential expression
    if ncells.shape[1] < 2:
        __import__("ipdb").set_trace()
        return pd.DataFrame()

    # You can only compute differential expression if both disease and control are present
    ncells = ncells.loc[ncells.all(axis=1)]
    comparisons = ncells.index
    comparisons = pd.DataFrame.from_records(comparisons, columns=comparisons.names)

    result = []
    # Group by dataset_id for speed
    for dataset_id, obs in comparisons.groupby("dataset_id"):
        # All the entries in obs are guaranteed to have nonzero cells for both normal and disease
        approx = scquill.Approximation.read_h5(get_dataset_path(dataset_id))
        # NOTE:: we should binarize consistently
        adata = approx.to_anndata(
            groupby=["tissue_general", "cell_type", "disease"],
        )

        if feature is not None:
            adata = adata[:, feature]

        obs_names = obs[["tissue_general", "cell_type"]].agg("\t".join, axis=1).values

        # Split between normal and disease
        adata1 = adata[adata.obs["disease"] == "normal"]
        adata1.obs.index = pd.Index(
            adata1.obs[["tissue_general", "cell_type"]].agg("\t".join, axis=1).values,
            name="category",
        )

        adata_notnormal = adata[adata.obs["disease"] != "normal"]
        disease_states = adata_notnormal.obs["disease"].unique()
        for disease in disease_states:
            adata2 = adata_notnormal[adata_notnormal.obs["disease"] == disease]
            adata2.obs.index = pd.Index(
                adata2.obs[["tissue_general", "cell_type"]]
                .agg("\t".join, axis=1)
                .values,
                name="category",
            )
            # Ok, now there might be cell types that are missing in one of the conditions, so we need to
            # use the obs variable to filter these adata
            obs_namesd = list(set(obs_names) & set(adata2.obs_names))
            adata1d = adata1[obs_namesd]
            adata2 = adata2[obs_namesd]

            result_dataset_id = _diff_exp_adatas(
                adata1d,
                adata2,
                number,
                dataset_id,
                method=method,
                feature=feature,
            )
            result.extend(result_dataset_id)

    result = pd.DataFrame(result).sort_values("metric", ascending=False)

    return result


def _diff_exp_adatas(
    adata1d, adata2, number, dataset_id, method="delta_fraction", feature=None
):
    """Compute differential expression between two approximated AnnData objects."""
    result = []
    if method == "delta_fraction":
        frac1 = adata1d.layers["fraction"]
        frac2 = adata2.layers["fraction"]
        delta_fraction = frac2 - frac1

        if feature is None:
            idx_bot = np.argpartition(delta_fraction.ravel(), number)[:number]
            idx_bot = np.unravel_index(idx_bot, delta_fraction.shape)
            idx_top = np.argpartition(delta_fraction.ravel(), -number)[-number:]
            idx_top = np.unravel_index(idx_top, delta_fraction.shape)
            idx_dict = {"down": idx_bot, "up": idx_top}
        else:
            idx = (
                np.arange(delta_fraction.shape[0]),
                np.zeros(delta_fraction.shape[0], dtype=int),
            )
            idx_dict = {"any": idx}

        for idxtype, idx in idx_dict.items():
            for i, j in zip(*idx):
                row1 = adata1d.obs.iloc[i]
                row2 = adata2.obs.iloc[i]
                res = {
                    "dataset_id": dataset_id,
                    "tissue_general": row1["tissue_general"],
                    "cell_type": row1["cell_type"],
                    "comparison": "disease vs. normal",
                    "condition": row2["disease"],
                    "condition_baseline": row1["disease"],
                    "regulation": idxtype,
                    "gene": adata1d.var_names[j],
                    "unit": "cptt",
                    "normal_expr": adata1d.layers["average"][i, j],
                    "disease_expr": adata2.layers["average"][i, j],
                    "normal_fraction": adata1d.layers["fraction"][i, j],
                    "disease_fraction": adata2.layers["fraction"][i, j],
                    "metric": delta_fraction[i, j],
                }
                result.append(res)
    else:
        raise NotImplementedError
    return result
