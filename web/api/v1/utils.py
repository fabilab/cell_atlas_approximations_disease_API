"""Utility functions for the API."""
import pandas as pd
from models.metadata import get_metadata
from models.exceptions import (
    NoContrastingConditionsInADatasetError,
    ParamsConflictError,
    UniqueIdNotFoundError,
)

# FIXME: there is a logical fallacy in allowing both unique_ids and metadta filters, deal with it at some point
def get_filter_kwargs(args, columns):
    """Return a dictionary of obs metadata keyword arguments to be used as database filters."""
    kwargs = {}
    
    unique_ids = args.get("unique_ids")
    
    if unique_ids:
        # Allow 'features' but not other metadata filters for specific functions, for example: average, dotplot ...
        allowed_keys = {"features", "unique_ids"}
        extra_keys = set(args.keys()) - allowed_keys
        if extra_keys:
            raise ParamsConflictError(
                msg=f"You can specify either unique_ids or metadata filters, not both"
            )
       
        unique_ids_list = unique_ids.split(",")
        
        # retrieve metadata for the given unique ids
        metadata = get_metadata()
        metadata_rows = metadata[metadata["unique_id"].isin(unique_ids_list)]

        if metadata_rows.empty:
            raise UniqueIdNotFoundError(
                msg=f"None of the provided unique IDs were found in metadata",
                unique_ids=unique_ids_list 
            )
        
        # Ensure each unique_id has a matching row with the same attributes except for disease.
        # This check is required for differential analysis, which needs both a disease and a normal (baseline) condition.
        for _, row in metadata_rows.iterrows():
            matching_rows = metadata[
                (metadata["dataset_id"] == row["dataset_id"]) &
                (metadata["cell_type"] == row["cell_type"]) &
                (metadata["tissue_general"] == row["tissue_general"]) &
                (metadata["development_stage_general"] == row["development_stage_general"]) &
                (metadata["sex"] == row["sex"]) &
                # Ensure one of them is "normal", in some edge cases, a dataset contains more than one disease conditions bit without normal, e.g: uids=8a859cb4f410fe3cb5a482896522b97b
                ((metadata["disease"] != row["disease"]) &
                ((metadata["disease"] == "normal") | (row["disease"] == "normal"))) 
            ]

            if matching_rows.empty:
                raise NoContrastingConditionsInADatasetError(
                    msg=(
                        "Differential analysis cannot be performed."
                        f"No contrasting condition found for unique_id '{row['unique_id']}' in dataset '{row['dataset_id']}'. "
                    ),
                    filters={"provided_unique_ids": unique_ids_list}
                )
        
        # If unique_ids include only normal samples, find corresponding disease cases
        normal_rows = metadata_rows[metadata_rows["disease"] == "normal"]
        if not normal_rows.empty:
            # Pair normal samples with disease samples from the same category
            disease_rows = metadata.merge(
                normal_rows,
                on=["dataset_id", "cell_type", "tissue_general", "development_stage_general", "sex"],
                suffixes=("", "normal"),
            ).query("disease != 'normal'")[metadata.columns]
            
            metadata_rows_filtered = metadata_rows[metadata_rows["disease"] != "normal"]

            metadata_rows = pd.concat([metadata_rows_filtered, disease_rows], ignore_index=True)
            
        # Extract metadata filters
        for column in columns:
            if column == 'tissue':
                column = 'tissue_general'
            
            if column == 'development_stage':
                column = 'development_stage_general'
            if column in metadata_rows.columns:
                kwargs[column] = metadata_rows[column].unique().tolist()

        return _clean_metadata_kwargs(kwargs)  
    
    # Handle non-unique_ids filtering
    for column in columns:
        value = args.get(column, default="", type=str)
        if value != "":
            # kwargs[column] = value
            # Convert comma-separated values into a list for multi-value filters
            kwargs[column] = value.split(",") if "," in value else value

    return _clean_metadata_kwargs(kwargs) 


def _clean_metadata_kwargs(kwargs):

    """Clean metadata filters"""
    if kwargs is None:
        return {}
    
    if "sex" in kwargs:
        value = _normalise_sex_string(kwargs.pop("sex"))
        if value is not None:
            kwargs["sex"] = value

    # NOTE: This is a temporary fix re renaming of column
    if "development_stage" in kwargs:
        kwargs["development_stage_general"] = kwargs.pop("development_stage")
    
    if "tissue" in kwargs:
        kwargs["tissue_general"] = kwargs.pop("tissue")
        
    if "unique_ids" in kwargs:
        if len(kwargs) > 1:
            raise ValueError(
            "You can only specify either unique_ids or metadata filters, not both"
        )
        unique_ids_str = kwargs.pop("unique_ids")
        # Skip if empty string
        if unique_ids_str:
            kwargs["unique_ids"] = unique_ids_str.split(",")

    return kwargs


def _normalise_sex_string(value):
    # Handle case where value is a list
    if isinstance(value, list):
        return [_normalise_sex_string(v) for v in value]  # Recursively normalize each value

    not_prefix, value = value.startswith("!"), value.lstrip("!")
    if value.lower().startswith("f"):
        value = "female"
    elif value.lower().startswith("m"):
        value = "male"
    elif value.lower() in ["any", "all"]:
        value = None
    if value is not None and not_prefix:
        value = f"!{value}"
    return value


def get_groupby_args(args, default=None):
    """Return a list of groupby variables."""
    groupby = args.get("groupby", None, type=str)
    if groupby is not None:
        groupby = groupby.replace(" ", "").split(",")
    else:
        groupby = default
    return groupby


def clean_feature_string(featurestring):
    """Clean feature string."""
    return featurestring.replace(" ", "").split(",")
