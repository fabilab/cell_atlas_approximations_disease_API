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
        keyword (str): Keyword to filter diseases.
        dataset_id (str): Identifier for the dataset.

    Returns:
        result (list): List of dictionaries containing the computed results.
"""
def compute_diff_cell_abundance(adata, keyword, dataset_id):
    df_obs = adata.obs
    filtered_obs = df_obs[df_obs['disease'].str.contains(keyword, case=False) | (df_obs['disease'] == 'normal')]
    disease_name = filtered_obs[filtered_obs['disease'].str.contains(keyword, case=False)]['disease'].unique()[0]
    result = []
    total_split_counts = filtered_obs.groupby('disease')['cell_count'].sum().to_dict()
    
    for cell_type in filtered_obs.cell_type.unique():
        normal_count = filtered_obs[(filtered_obs["cell_type"] == cell_type) & (filtered_obs["disease"] == 'normal')]["cell_count"].sum()
        disease_count = filtered_obs[(filtered_obs["cell_type"] == cell_type) & (filtered_obs["disease"] != 'normal')]["cell_count"].sum()
        total_count = normal_count + disease_count
        normal_pct = (normal_count / total_count) * 100
        disease_pct = (disease_count / total_count) * 100
        delta_fraction = disease_pct - normal_pct

        # OrderedDict is used to ensure the dictionary is shown in the same order as defined
        result.append(OrderedDict([
            ("disease", disease_name),
            ("dataset_id", dataset_id),
            ("cell_type", cell_type),
            ("normal_count", normal_count),
            ("disease_count", disease_count),
            ("total_count", total_count),
            ("normal_pct", normal_pct),
            ("disease_pct", disease_pct),
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

def compute_diff_expression(adata, keyword, dataset_id, N, cell_type_keyword):
    expression_data = adata.layers['average']
    fraction_data = adata.layers['fraction']
    filtered_obs = adata.obs[adata.obs['disease'].str.contains(keyword, case=False) | (adata.obs['disease'] == 'normal')]
    disease_name = filtered_obs[filtered_obs['disease'].str.contains(keyword, case=False)]['disease'].unique()[0]  # Extract the disease name
    cell_types = filtered_obs['cell_type'].unique()
    disease_status = filtered_obs['disease'].unique()

    result = []
    for cell_type in cell_types:
        if cell_type_keyword.lower() not in cell_type.lower():
            continue

        for status in disease_status:
            if status == 'normal':
                continue

            normal_idx = np.where((filtered_obs['cell_type'] == cell_type) & (filtered_obs['disease'] == 'normal'))[0][0]
            disease_idx = np.where((filtered_obs['cell_type'] == cell_type) & (filtered_obs['disease'] == status))[0][0]

            expr_normal = expression_data[normal_idx, :]
            expr_disease = expression_data[disease_idx, :]
            frac_normal = fraction_data[normal_idx, :]
            frac_disease = fraction_data[disease_idx, :]

            log2_fc = np.log2((expr_disease + 1) / (expr_normal + 1))
            
            if N >= log2_fc.shape[0]:
                N = log2_fc.shape[0] - 1

            top_up_indices = np.argsort(log2_fc)[-N:][::-1]
            top_down_indices = np.argsort(log2_fc)[:N]

            for idx in top_up_indices.tolist() + top_down_indices.tolist():
                regulation = 'up' if log2_fc[idx] > 0 else 'down'
                result.append(OrderedDict([
                    ("disease", disease_name),
                    ("dataset_id", dataset_id),
                    ("condition", status),
                    ("cell_type", cell_type),
                    ("regulation", regulation),
                    ("feature_name", adata.var_names[idx]),
                    ("expr_normal", expr_normal[idx]),
                    ("expr_disease", expr_disease[idx]),
                    ("frac_normal", frac_normal[idx]),
                    ("frac_disease", frac_disease[idx]),
                    ("log2_fc", log2_fc[idx]),
                    ("delta_change", (expr_disease[idx] - expr_normal[idx]))
                ]))

    return convert_to_python_types(result)