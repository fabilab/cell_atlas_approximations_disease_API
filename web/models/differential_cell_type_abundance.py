import pandas as pd

from models.metadata import get_metadata


def get_diff_cell_abundance(**filters):
    """
    Computes the differential cell type abundance.

    Parameters:
        **filters: Arbitrary keyword arguments to filter the metadata.

    Returns:
        result (list): List of dictionaries containing the computed results.
    """
    meta = get_metadata(**filters)

    result = []
    for dataset_id, obs in meta.groupby("dataset_id"):
        # Check if both disease and control are present in the dataset
        disease_states = obs["disease"].unique()
        if "normal" not in disease_states or len(disease_states) < 2:
            continue

        table = (
            obs.groupby(["cell_type", "disease"])["cell_count"]
            .sum()
            .unstack("disease", fill_value=0)
        )
        table.rename(columns={"normal": "ncell_normal"}, inplace=True)
        table["frac_normal"] = table["ncell_normal"] / table["ncell_normal"].sum()
        table["dataset_id"] = dataset_id
        for state in disease_states:
            if state == "normal":
                continue
            # Slight optimisation whenever there is only one disease state
            # avoids a copy
            if len(disease_states) == 2:
                table_state = table
            else:
                table_state = table[
                    ["dataset_id", state, "ncell_normal", "frac_normal"]
                ].copy()
            table_state.rename(columns={state: "ncell_disease"}, inplace=True)
            table_state["frac_disease"] = (
                table_state["ncell_disease"] / table_state["ncell_disease"].sum()
            )
            table_state["delta_frac"] = (
                table_state["frac_disease"] - table_state["frac_normal"]
            )
        table_state["disease"] = state
        result.append(table_state)

    result = pd.concat(result)

    return result
