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
    
    # Set numerical indices on the original adata.obs
    adata.obs['numerical_index'] = np.arange(adata.obs.shape[0])
    
    # Filter the obs to only include rows with the disease_keyword or 'normal'
    filtered_obs = adata.obs[adata.obs['disease'].str.contains(disease_keyword, case=False) | (adata.obs['disease'] == 'normal')]
    disease_name = filtered_obs[filtered_obs['disease'].str.contains(disease_keyword, case=False)]['disease'].unique()[0]
    
    # Convert filtered_obs['numerical_index'] to a numpy array of integer indices
    integer_indices = filtered_obs['numerical_index'].to_numpy()

    # Filter the expression and fraction data using the filtered observations' integer indices
    filtered_expression_data = expression_data[integer_indices, :]
    filtered_fraction_data = fraction_data[integer_indices, :]
    
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

            # Get the numerical indices for normal and disease
            normal_idx = filtered_obs[(filtered_obs['cell_type'] == cell_type) & (filtered_obs['disease'] == 'normal')]['numerical_index'].values[0]
            disease_idx = filtered_obs[(filtered_obs['cell_type'] == cell_type) & (filtered_obs['disease'] == status)]['numerical_index'].values[0]

            # Ensure normal_idx and disease_idx are within bounds
            if normal_idx >= len(filtered_expression_data) or disease_idx >= len(filtered_expression_data):
                continue

            normal_expr = filtered_expression_data[normal_idx, :]
            disease_expr = filtered_expression_data[disease_idx, :]
            normal_fraction = filtered_fraction_data[normal_idx, :]
            disease_fraction = filtered_fraction_data[disease_idx, :]
            delta_fraction = disease_fraction - normal_fraction

            log2_fc = np.log2((disease_expr + 1) / (normal_expr + 1))

            # ensure top_n does not exceed the number of features
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
                    ("log_transformed", log_transformed),
                    ("log2_fc", log2_fc[idx]),
                    ("normal_fraction", normal_fraction[idx]),
                    ("disease_fraction", disease_fraction[idx]),
                    ("delta_fraction", delta_fraction[idx])
                ]))

    return convert_to_python_types(result)


def get_metadata(disease_keyword, metadata):
    
    result = []
    for disease, unique_id in zip(metadata['diseases'], metadata['unique_ids']):
        if disease_keyword.lower() in disease.decode('utf-8').lower():
            result.append(OrderedDict([
                ('unique_id', unique_id.decode('utf-8')),
                ('disease', disease.decode('utf-8')),
                ('cell_type_number', len(metadata['cell_types'])),
                ('dataset_id', metadata['dataset_id']),
                ('unit', metadata['unit'].decode('utf-8') if isinstance(metadata['unit'], bytes) else metadata['unit']),
                ('log_transformed', metadata['log_transformed']),
                ('has_normal_baseline', metadata['has_normal_baseline'])
            ]))

    return convert_to_python_types(result)