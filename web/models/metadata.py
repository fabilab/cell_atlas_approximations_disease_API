import re
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
        # Boolean AND between filters
        keep &= metadata[key].str.contains(value, case=False)
    return metadata.loc[keep]
