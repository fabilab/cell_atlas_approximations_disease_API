import re
import numpy as np
import pandas as pd
import hashlib

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
    """Gnerate SHA256 hash for a given row"""
    row_string = ",".join(map(str, row))  # Convert row values to a single string
    m = hashlib.md5()
    m.update(row_string.encode("utf-8"))  # Update the hash with encoded row string
    return m.hexdigest()  # Generate MD5 hash

def load_metadata():
    """Cache metadata from the manifest file."""
    global metadata
    obs = pd.read_csv(config["paths"]["obs_metadata_file"])
    metadata = obs
    
    # Add unique_id column to metadata. Generate unique_id ONCE when starting the server and store it
    # If already stored, do not touch
    if "unique_id" not in obs.columns:
        obs["unique_id"] = obs.apply(generate_unique_id, axis=1)


def get_metadata(**filters):
    """Get metadata that filfill all given filters."""
    # Lazy loading of metadata cache upon first call
    if metadata is None:
        load_metadata()

    keep = pd.Series(np.ones(len(metadata), dtype=bool), index=metadata.index)
    for key, value in filters.items():
        if key == "unique_ids":
            continue
        
        # Boolean OR
        invert, value = value.startswith("!"), value.lstrip("!")

        # Boolean AND between filters
        if key == "sex":
            # "male" is a substring of "female"
            keep_key = metadata[key] == value
        elif key == "cell_type":
            value = value.strip()
            # If the cell_type contains a space (e.g., 't cell'), use strict matching to ensure it's a standalone term
            # e.g., 'goblet cell' and 'mast cell' shouldn't be included as a type of 't cell'
            if " " in value:
                value = rf"(?i)(?<!\w){value}(?!\w)"
            # If the cell_type is a single word (e.g., 'neuron'), match as a suffix to include subtypes like 'interneuron'
            # but ensure the word is not preceded by letters to avoid false positives
            else:
                value = rf"(?i)(?<!\w)(?:\w*{value})(?!\w)"
            keep_key = metadata[key].str.contains(value, na=False, regex=True)
        else:
            keep_key = metadata[key].str.contains(value, case=False)

        if invert:
            keep_key = ~keep_key
        keep &= keep_key

    filtered_metadata = metadata.loc[keep]

    # Handle errors separately for each parameter
    if "disease" in filters and metadata["disease"].str.contains(filters["disease"], case=False).sum() == 0:
        raise DiseaseNotFoundError(
            msg=f"No disease found that matches '{filters['disease']}'.",
            disease=filters["disease"]
        )

    if "cell_type" in filters:
        value = filters["cell_type"].strip()
        if " " in value:
            value = rf"(?i)(?<!\w){value}(?!\w)"
        else:
            value = rf"(?i)(?<!\w)(?:\w*{value})(?!\w)"
        if metadata["cell_type"].str.contains(value, na=False, regex=True).sum() == 0:
            raise CellTypeNotFoundError(
                msg=f"No cell type found that matches '{filters['cell_type']}'.",
                cell_type=filters["cell_type"]
            )
        
    if "tissue" in filters and metadata["tissue"].str.contains(filters["tissue"], case=False).sum() == 0:
        raise TissueNotFoundError(
            msg=f"No tissue found that matches '{filters['tissue']}'.",
            tissue=filters["tissue"]
        )

    if "development_stage_general" in filters and metadata["development_stage_general"].str.contains(filters["development_stage_general"], case=False).sum() == 0:
        raise DevelopmentStageNotFoundError(
            msg=f"No development stage found that matches '{filters['development_stage_general']}'.",
            development_stage_general=filters["development_stage_general"]
        )
    
    if filtered_metadata.empty:
        raise NoMatchingDatasetsError(
            msg="No datasets found that satisfy the requested filters.",
            filters=filters
        )

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
