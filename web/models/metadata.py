import re
import numpy as np
import pandas as pd
from config import configuration as config
from models.utils import (
    convert_to_python_types,
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

    return metadata.loc[keep]


def get_metadata_with_normal(**filters):
    """Get metadata that fulfill all given filters, including normal samples."""
    if "disease" not in filters or filters["disease"] == "normal":
        return get_metadata(**filters)

    # Disease was among the filters, so we need to get both disease and normal samples
    meta_disease = get_metadata(**filters)

    filters_normal = filters.copy()
    filters_normal["disease"] = "normal"
    meta_normal = get_metadata(**filters_normal)

    meta_joint = pd.concat([meta_disease, meta_normal])
    return meta_joint
