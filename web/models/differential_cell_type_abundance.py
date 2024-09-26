import pandas as pd

from models.metadata import get_metadata


def get_diff_cell_abundance(**filters):
    """
    Computes the differential cell type abundance for a given dataset.

    Parameters:
        adata (AnnData): data obtained from scquill
        disease_keyword (str): Keyword to filter diseases.
        dataset_id (str): Identifier for the dataset.

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

        # NOTE: this should be improved
        meta["disease_bin"] = meta["disease"] != "normal"

        table = (
            meta.groupby(["cell_type", "disease_bin"])["cell_count"]
            .sum()
            .unstack("disease_bin", fill_value=0)
        )
        table.rename(
            columns={True: "ncells_disease", False: "ncells_normal"}, inplace=True
        )
        table["frac_disease"] = table["ncells_disease"] / table["ncells_disease"].sum()
        table["frac_normal"] = table["ncells_normal"] / table["ncells_normal"].sum()
        table["delta_frac"] = table["frac_disease"] - table["frac_normal"]
        table["dataset_id"] = dataset_id
        table["ncells_total_disease"] = table["ncells_disease"].sum()
        table["ncells_total_normal"] = table["ncells_normal"].sum()
        result.append(table)

    result = pd.concat(result)

    return result

