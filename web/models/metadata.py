import numpy as np
import pandas as pd

from config import configuration as config
from models.baseline import get_differential_baseline

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

    # Map API parameters to database columns
    if 'tissue' in filters:
        filters['tissue_general'] = filters.pop('tissue')

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

    return metadata.loc[keep]


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
