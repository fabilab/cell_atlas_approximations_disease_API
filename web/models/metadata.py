import numpy as np
import pandas as pd
import hashlib
import time

from config import configuration as config
from models.baseline import get_differential_baseline
from models.exceptions import (
    DiseaseNotFoundError,
    CellTypeNotFoundError,
    TissueNotFoundError,
    DevelopmentStageNotFoundError,
    NoMatchingDatasetsError,
)

metadata = None

def generate_unique_id(row):
    """Generate md5 hash (faster than 256) for a given row"""
    row_string = ",".join(map(str, row))  # Convert row values to a single string
    m = hashlib.md5()
    m.update(row_string.encode("utf-8"))  # Update the hash with encoded row string
    return m.hexdigest()  # Generate MD5 hash

def load_metadata():
    """Cache metadata from the manifest file."""
    global metadata
    obs = pd.read_csv(config["paths"]["obs_metadata_file"])
    metadata = obs
    
    # Measure time for generating unique IDs
    unique_id_start_time = time.time()
    # Add unique_id column to metadata. Generate unique_id ONCE and store it
    obs["unique_id"] = obs.apply(generate_unique_id, axis=1)
    unique_id_end_time = time.time()
    print(f"Time taken to generate unique IDs: {unique_id_end_time - unique_id_start_time:.4f} seconds")

def get_metadata(**filters):
    """Get metadata that filfill all given filters."""
    # Lazy loading of metadata cache upon first call
    if metadata is None:
        load_metadata()

    keep = pd.Series(np.ones(len(metadata), dtype=bool), index=metadata.index)
    for key, value in filters.items():
        # Boolean OR
        if isinstance(value, str):
            invert, value = value.startswith("!"), [value.lstrip("!")]

        else:
            invert = False
            value = [v.lstrip("!") for v in value]

        # Boolean AND between filters
        if key == "sex":
            # If value is a single string, use `==`
            if isinstance(value, str) or len(value) == 1:
                keep_key = metadata[key] == value[0] if isinstance(value, list) else metadata[key] == value
            else:
                # If value is a list, use `.isin()` to prevent shape mismatch
                keep_key = metadata[key].isin(value)
    
        else:
            # e.g: ["B cell", "T cell"]
            if isinstance(value, list):
                keep_key = pd.Series(np.zeros(len(metadata), dtype=bool), index=metadata.index)
                for v in value:
                    keep_key |= metadata[key].str.contains(v.strip(), case=False)
            else:
                keep_key = metadata[key].str.contains(value, case=False)

        if invert:
            keep_key = ~keep_key
        keep &= keep_key

    filtered_metadata = metadata.loc[keep]

    # Handle errors separately for each parameter
    if "disease" in filters:
        if (isinstance(filters["disease"], str) and metadata["disease"].str.contains(filters["disease"], case=False).sum() == 0) \
        or (isinstance(filters["disease"], list) and all(metadata["disease"].str.contains(v, case=False).sum() == 0 for v in filters["disease"])):
            raise DiseaseNotFoundError(
                msg=f"No disease found that matches '{filters['disease']}'.",
                disease=filters["disease"]
            )

    if "cell_type" in filters:
        if (isinstance(filters["cell_type"], str) and metadata["cell_type"].str.contains(filters["cell_type"], case=False).sum() == 0) \
        or (isinstance(filters["cell_type"], list) and all(metadata["cell_type"].str.contains(v, case=False).sum() == 0 for v in filters["cell_type"])):
            raise CellTypeNotFoundError(
                msg=f"No cell type found that matches '{filters['cell_type']}'.",
                cell_type=filters["cell_type"]
            )
        
    if "tissue" in filters:
        if (isinstance(filters["tissue"], str) and metadata["tissue"].str.contains(filters["tissue"], case=False).sum() == 0) \
        or (isinstance(filters["tissue"], list) and all(metadata["tissue"].str.contains(v, case=False).sum() == 0 for v in filters["tissue"])):
            raise TissueNotFoundError(
                msg=f"No tissue found that matches '{filters['tissue']}'.",
                tissue=filters["tissue"]
            )

    if "development_stage_general" in filters:
        if (isinstance(filters["development_stage_general"], str) and metadata["development_stage_general"].str.contains(filters["development_stage_general"], case=False).sum() == 0) \
        or (isinstance(filters["development_stage_general"], list) and all(metadata["development_stage_general"].str.contains(v, case=False).sum() == 0 for v in filters["development_stage_general"])):
            raise DevelopmentStageNotFoundError(
                msg=f"No development stage found that matches '{filters['development_stage_general']}'.",
                development_stage_general=filters["development_stage_general"]
            )
    
    if filtered_metadata.empty:
        raise NoMatchingDatasetsError(
            msg="No datasets found that satisfy the requested filters.",
            filters=filters
        )
        
    # Generate unique_id AFTER filtering
    # filtered_metadata["unique_id"] = filtered_metadata.apply(generate_unique_id, axis=1)

    # Move `unique_id` to the first column
    cols = ["unique_id"] + [col for col in filtered_metadata.columns if col != "unique_id"]
    filtered_metadata = filtered_metadata[cols]

    return filtered_metadata


def get_metadata_with_baseline(differential_axis, **filters):
    """Get metadata that fulfill all given filters, including baseline samples."""
    baseline = get_differential_baseline(differential_axis)
    if differential_axis not in filters or filters[differential_axis] == baseline:
        return get_metadata(**filters)

    # The differential was among the filters, so we need to get both state and baseline samples
    meta_state = get_metadata(**filters)

    filters_baseline = filters.copy()
    filters_baseline[differential_axis] = baseline
    meta_baseline = get_metadata(**filters_baseline)

    meta_joint = pd.concat([meta_state, meta_baseline])
    return meta_joint
