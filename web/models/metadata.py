import numpy as np
import pandas as pd

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


def load_metadata():
    """Cache metadata from the manifest file."""
    global metadata
    obs = pd.read_csv(config["paths"]["obs_metadata_file"])
    metadata = obs


def get_metadata(**filters):
    """Get metadata that filfill all given filters."""
    # Lazy loading of metadata cache upon first call
    if metadata is None:
        load_metadata()

    keep = pd.Series(np.ones(len(metadata), dtype=bool), index=metadata.index)
    for key, value in filters.items():
        # Boolean OR
        invert, value = value.startswith("!"), value.lstrip("!")

        # Boolean AND between filters
        if key == "sex":
            # "male" is a substring of "female"
            keep_key = metadata[key] == value
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

    if "cell_type" in filters and metadata["cell_type"].str.contains(filters["cell_type"], case=False).sum() == 0:
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
