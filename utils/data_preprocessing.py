# utils/data_processing.py

import numpy as np
import pandas as pd
from collections import OrderedDict

"""
result (list): List of dictionaries containing results.
Returns (list): List of dictionaries with values converted to Python types.
This function is used to fix the error of: TypeError: Object of type int64 is not JSON serializable
"""
def convert_to_python_types(result):
    for entry in result:
        for key, value in entry.items():
            if isinstance(value, (np.integer, np.int64)):
                entry[key] = int(value)
            elif isinstance(value, (np.floating, np.float64)):
                entry[key] = float(value)
            elif isinstance(value, (np.ndarray, np.generic)):
                entry[key] = value.tolist()
    return result


"""
    Computes the differential cell type abundance for a given dataset.

    Parameters:
        adata (AnnData): data obtained from scquill
        disease_keyword (str): Keyword to filter diseases.
        dataset_id (str): Identifier for the dataset.

    Returns:
        result (list): List of dictionaries containing the computed results.
"""
def compute_diff_cell_abundance(adata, disease_keyword, dataset_id):
    df_obs = adata.obs
    filtered_obs = df_obs[df_obs['disease'].str.contains(disease_keyword, case=False) | (df_obs['disease'] == 'normal')]
    disease_name = filtered_obs[filtered_obs['disease'].str.contains(disease_keyword, case=False)]['disease'].unique()[0]
    result = []
    
    for cell_type in filtered_obs.cell_type.unique():
        normal_count = filtered_obs[(filtered_obs["cell_type"] == cell_type) & (filtered_obs["disease"] == 'normal')]["cell_count"].sum()
        disease_count = filtered_obs[(filtered_obs["cell_type"] == cell_type) & (filtered_obs["disease"] != 'normal')]["cell_count"].sum()
        total_count = normal_count + disease_count
        normal_fraction = normal_count / total_count
        disease_fraction = disease_count / total_count
        delta_fraction = disease_fraction - normal_fraction

        # OrderedDict is used to ensure the dictionary is shown in the same order as defined
        result.append(OrderedDict([
            ("disease", disease_name),
            ("dataset_id", dataset_id),
            ("cell_type", cell_type),
            ("comparison", "disease vs. normal"),
            ("condition", "not normal"),
            ("condition_baseline", "normal"),
            ("normal_count", normal_count),
            ("disease_count", disease_count),
            ("total_count", total_count),
            ("normal_fraction", normal_fraction),
            ("disease_fraction", disease_fraction),
            ("delta_fraction", delta_fraction)
        ]))
    
    return convert_to_python_types(result)


"""
    Computes the differential cell type abundance for a given dataset.

    Parameters:
        adata (AnnData): data obtained from scquill
        keyword (str): Keyword to filter diseases.
        dataset_id (str): Identifier for the dataset.

    Returns:
        result (list): List of dictionaries containing the computed results.
"""

def compute_diff_expression(adata, disease_keyword, dataset_id, metadata, N, cell_type_keyword):
    expression_data = adata.layers['average']
    fraction_data = adata.layers['fraction']
    filtered_obs = adata.obs[adata.obs['disease'].str.contains(disease_keyword, case=False) | (adata.obs['disease'] == 'normal')]
    disease_name = filtered_obs[filtered_obs['disease'].str.contains(disease_keyword, case=False)]['disease'].unique()[0]  # Extract the disease name
    cell_types = filtered_obs['cell_type'].unique()
    disease_status = filtered_obs['disease'].unique()
    unit = metadata['unit']
    log_transformed = metadata['log_transformed']
    
    result = []
    for cell_type in cell_types:
        if cell_type_keyword.lower() not in cell_type.lower():
            continue

        for status in disease_status:
            if status == 'normal':
                continue

            normal_idx = np.where((filtered_obs['cell_type'] == cell_type) & (filtered_obs['disease'] == 'normal'))[0][0]
            disease_idx = np.where((filtered_obs['cell_type'] == cell_type) & (filtered_obs['disease'] == status))[0][0]

            normal_expr = expression_data[normal_idx, :]
            disease_expr = expression_data[disease_idx, :]
            normal_fraction = fraction_data[normal_idx, :]
            disease_fraction = fraction_data[disease_idx, :]
            delta_fraction = disease_fraction - normal_fraction

            log2_fc = np.log2((disease_expr + 1) / (normal_expr + 1))

            if N >= log2_fc.shape[0]:
                N = log2_fc.shape[0] - 1

            # Sort by delta fraction
            top_up_indices = np.argsort(delta_fraction)[-N:][::-1]
            top_down_indices = np.argsort(delta_fraction)[:N]

            # Combine up and down indices with regulation labels
            combined_indices = [(idx, 'up') for idx in top_up_indices] + [(idx, 'down') for idx in top_down_indices]

            for idx, regulation in combined_indices:
                result.append(OrderedDict([
                    ("disease", disease_name),
                    ("dataset_id", dataset_id),
                    ("cell_type", cell_type),
                    ("comparison", "disease vs. normal"),
                    ("condition", "not normal"),
                    ("condition_baseline", "normal"),
                    ("regulation", regulation),
                    ("feature_name", adata.var_names[idx]),
                    ("unit", unit),
                    ("normal_expr", normal_expr[idx]),
                    ("disease_expr", disease_expr[idx]),
                    ("log2_fc", log2_fc[idx]),
                    ("normal_fraction", normal_fraction[idx]),
                    ("disease_fraction", disease_fraction[idx]),
                    ("delta_fraction", delta_fraction[idx])
                ]))

    return convert_to_python_types(result)