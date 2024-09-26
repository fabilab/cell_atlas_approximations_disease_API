"""Module for handling paths."""

import pathlib
from config import configuration as config


def get_dataset_path(dataset_id):
    """Get the path to the dataset file."""
    atlas_path = pathlib.Path(config["paths"]["compressed_atlas"])
    return atlas_path / f"{dataset_id}.h5"
