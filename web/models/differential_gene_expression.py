import numpy as np
import pandas as pd
import scquill

from models.paths import get_dataset_path
from models.baseline import get_differential_baseline
from models.metadata import get_metadata_with_baseline
from models.exceptions import (
    NoContrastingConditionsInADatasetError
)


_differential_baselines = {
    "disease": "normal",
    "sex": "female",
    "age": "adult",
}


def get_diff_expression(
    differential_axis="disease",
    groupby=None,
    number=10,
    feature=None,
    method="delta_fraction",
    **filters,
):
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

    if groupby is None:
        groupby = []
    # Ensure the comparison group has the same cell type and tissue (by default)
    if "cell_type" not in groupby:
        groupby = ["cell_type"] + groupby
        for i, name in enumerate(groupby):
            if name == "tissue":
                groupby[i] = "tissue_general"

    if differential_axis in groupby:
        raise ValueError(
            f"{differential_axis} cannot be a groupby variable and the differential axis at the same time"
        )

    baseline = get_differential_baseline(differential_axis)
    meta = get_metadata_with_baseline(
        differential_axis,
        **filters,
    )
    if differential_axis not in meta.columns:
        raise ValueError(f"Metadata does not contain {differential_axis}")

    if len(meta) == 0:
        return pd.DataFrame()

    result = []
    # NOTE: This grouping explicitely forbids cross-dataset comparisons
    # that's on purpose for now, re batch effects
    for dataset_id, obs in meta.groupby("dataset_id"):
        # Check if both focal group and baseline are present in the dataset
        differential_states = obs[differential_axis].unique()
        if baseline not in differential_states or len(differential_states) < 2:
            continue

        table = (
            obs.groupby(groupby + [differential_axis])["cell_count"]
            .sum()
            .unstack(differential_axis, fill_value=0)
        )
        cell_types_to_keep = table.index.get_level_values('cell_type').unique().tolist()
        table["dataset_id"] = dataset_id

        # All the entries in obs are guaranteed to have nonzero cells for both normal and disease
        approx = scquill.Approximation.read_h5(get_dataset_path(dataset_id))
        # NOTE:: we should binarize consistently
        adata = approx.to_anndata(
            groupby=groupby + [differential_axis],
        )
        
        # In the previous code, the calculation is perform without filtering cell type provided by the user,
        # There we  
        mask = adata.obs['cell_type'].isin(cell_types_to_keep)
    
        filtered_adata = adata[mask].copy()

        if feature is not None:
            filtered_adata = filtered_adata[:, feature]

        obs_names = obs[groupby].agg("\t".join, axis=1).values

        # Split between normal and disease
        adata_baseline = filtered_adata[filtered_adata.obs[differential_axis] == baseline]
        adata_baseline.obs.index = pd.Index(
            adata_baseline.obs[groupby].agg("\t".join, axis=1).values,
            name="category",
        )
        for state in differential_states:
            if state == baseline:
                continue
            adata_state = filtered_adata[filtered_adata.obs[differential_axis] == state]
            adata_state.obs.index = pd.Index(
                adata_state.obs[groupby].agg("\t".join, axis=1).values,
                name="category",
            )
            # Ok, now there might be cell types that are missing in one of the conditions, so we need to
            # use the obs variable to filter these adata
            # obs_namesd = list(set(obs_names) & set(adata_state.obs_names))
            
            obs_namesd = list(set(adata_baseline.obs_names) & set(adata_state.obs_names))
            adata1d = adata_baseline[obs_namesd]
            adata2 = adata_state[obs_namesd]

            result_dataset_id = _diff_exp_adatas(
                adata1d,
                adata2,
                number,
                method=method,
                feature=feature,
            )
            for res in result_dataset_id:
                res["dataset_id"] = dataset_id
                res["differential_axis"] = differential_axis
                res["state"] = state
                res["baseline"] = baseline
            result.extend(result_dataset_id)
    
    if len(result) == 0:
        raise NoContrastingConditionsInADatasetError(
            msg=f"No datasets can be found with more than one value along the '{differential_axis}' axis after filtering",
            filters=filters
        )
    
    result = pd.DataFrame(result).sort_values("metric", ascending=False)

    return result


def _diff_exp_adatas(adata1d, adata2, number, method="delta_fraction", feature=None):
    """Compute differential expression between two approximated AnnData objects."""
    result = []
    if method == "delta_fraction":
        frac1 = adata1d.layers["fraction"]
        frac2 = adata2.layers["fraction"]
        delta = frac2 - frac1
    elif method == "ratio_average":
        avg1 = adata1d.layers["average"]
        avg2 = adata2.layers["average"]
        delta = (avg2 + 1e-5) / (avg1 + 1e-5)
    else:
        raise NotImplementedError

    if feature is None:
        idx_bot = np.argpartition(delta.ravel(), number)[:number]
        idx_bot = np.unravel_index(idx_bot, delta.shape)
        idx_top = np.argpartition(delta.ravel(), -number)[-number:]
        idx_top = np.unravel_index(idx_top, delta.shape)
        idx_dict = {"down": idx_bot, "up": idx_top}
    else:
        idx = (
            np.arange(delta.shape[0]),
            np.zeros(delta.shape[0], dtype=int),
        )
        idx_dict = {"any": idx}

    for idxtype, idx in idx_dict.items():
        for i, j in zip(*idx):
            row1 = adata1d.obs.iloc[i]
            res = {
                "tissue_general": row1["tissue_general"],
                "cell_type": row1["cell_type"],
                "regulation": idxtype,
                "gene": adata1d.var_names[j],
                "unit": "cptt",
                "baseline_expr": adata1d.layers["average"][i, j],
                "state_expr": adata2.layers["average"][i, j],
                "baseline_fraction": adata1d.layers["fraction"][i, j],
                "state_fraction": adata2.layers["fraction"][i, j],
                "metric": delta[i, j],
            }
            result.append(res)
    return result
