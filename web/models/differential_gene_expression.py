import numpy as np
import pandas as pd
import scquill

from models.metadata import get_metadata
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

    meta = get_metadata(**filters)

    # NOTE: this should be improved
    meta["disease_bin"] = meta["disease"] != "normal"

    ncells = (
        meta.groupby(["dataset_id", "tissue_general", "cell_type", "disease_bin"])[
            "cell_count"
        ]
        .sum()
        .unstack("disease_bin", fill_value=0)
    )
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

        # Split between normal and disease
        adata1 = adata[adata.obs["disease"] == "normal"]
        adata2 = adata[adata.obs["disease"] != "normal"]

        # Align the two AnnData objects
        adata1.obs.index = pd.Index(
            adata1.obs[["tissue_general", "cell_type"]].agg("\t".join, axis=1).values,
            name="category",
        )
        # FIXME: this could lead to nonunique indices if there are multiple disease conditions
        adata2.obs.index = pd.Index(
            adata2.obs[["tissue_general", "cell_type"]].agg("\t".join, axis=1).values,
            name="category",
        )
        if adata2.obs_names.duplicated().any():
            continue
        obs_names = obs[["tissue_general", "cell_type"]].agg("\t".join, axis=1).values

        adata1 = adata1[obs_names]
        adata2 = adata2[obs_names]

        # Ok, now there might be cell types that are missing in one of the conditions, so we need to
        # use the obs variable to filter these adata

        if method == "delta_fraction":
            frac1 = adata1.layers["fraction"]
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
                    row1 = adata1.obs.iloc[i]
                    row2 = adata2.obs.iloc[i]
                    res = {
                        "dataset_id": dataset_id,
                        "tissue_general": row1["tissue_general"],
                        "cell_type": row1["cell_type"],
                        "comparison": "disease vs. normal",
                        "condition": row2["disease"],
                        "condition_baseline": row1["disease"],
                        "regulation": idxtype,
                        "gene": adata1.var_names[j],
                        "unit": "cptt",
                        "normal_expr": adata1.layers["average"][i, j],
                        "disease_expr": adata2.layers["average"][i, j],
                        "normal_fraction": adata1.layers["fraction"][i, j],
                        "disease_fraction": adata2.layers["fraction"][i, j],
                        "metric": delta_fraction[i, j],
                    }
                    result.append(res)
        else:
            raise NotImplementedError

    result = pd.DataFrame(result)

    return result
