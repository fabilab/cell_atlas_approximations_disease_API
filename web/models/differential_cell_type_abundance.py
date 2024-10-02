import pandas as pd

from models.metadata import get_metadata
from models.baseline import get_differential_baseline


def get_diff_cell_abundance(
    differential_axis="disease",
    groupby=None,
    **filters,
):
    """
    Computes the differential cell type abundance.

    Parameters:
        **filters: Arbitrary keyword arguments to filter the metadata.

    Returns:
        result (list): List of dictionaries containing the computed results.
    """
    if groupby is None:
        groupby = []
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
    meta = get_metadata(**filters)
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
        table["dataset_id"] = dataset_id
        table.rename(columns={baseline: "ncell_baseline"}, inplace=True)
        table["frac_baseline"] = table["ncell_baseline"] / table["ncell_baseline"].sum()
        for state in differential_states:
            if state == baseline:
                continue
            # Slight optimisation whenever there is only one disease state
            table_state = table[
                ["dataset_id", state, "ncell_baseline", "frac_baseline"]
            ]
            if len(differential_states) > 2:
                table_state = table_state.copy()
            table_state.rename(
                columns={state: f"ncell_{differential_axis}"}, inplace=True
            )
            table_state[f"frac_{differential_axis}"] = (
                table_state[f"ncell_{differential_axis}"]
                / table_state[f"ncell_{differential_axis}"].sum()
            )
            table_state["delta_frac"] = (
                table_state[f"frac_{differential_axis}"] - table_state["frac_baseline"]
            )
            table_state[differential_axis] = state
            table_state["baseline"] = baseline
            table_state.reset_index(inplace=True)
            result.append(table_state)

    result = pd.concat(result)

    # Reorder columns so that datset_id is the first column
    cols = result.columns.tolist()
    cols.remove("dataset_id")
    cols = ["dataset_id"] + cols
    result = result[cols]

    return result
