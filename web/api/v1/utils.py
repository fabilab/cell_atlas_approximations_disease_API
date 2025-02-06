"""Utility functions for the API."""
from models.metadata import get_metadata
from models.exceptions import (
    ParamsConflictError,
    UniqueIdNotFoundError,
)

# FIXME: there is a logical fallacy in allowing both unique_ids and metadta filters, deal with it at some point
def get_filter_kwargs(args, columns):
    """Return a dictionary of obs metadata keyword arguments to be used as database filters."""
    kwargs = {}
    
    # if unique ids are provided
    unique_ids = args.get("unique_ids")
    
    if unique_ids:
        if len(args) > 1:
            raise ParamsConflictError(
                msg=f"You can specify either unique_ids or metadata filters, not both"
            )
       
        unique_ids_list = unique_ids.split(",")
        
        # retrieve metadata for the given unique ids
        metadata = get_metadata()
        metadata_rows = metadata[metadata["unique_id"].isin(unique_ids_list)]
        # print(metadata_rows.columns)
        
        if metadata_rows.empty:
            raise UniqueIdNotFoundError(
                msg=f"None of the provided unique IDs were found in metadata",
                unique_ids=unique_ids_list 
            )
        

        # Extract fields from metadata rows and apply filtering
        for column in columns:
            if column in metadata_rows.columns:
                kwargs[column] = metadata_rows[column].unique().tolist()
         
        return _clean_metadata_kwargs(kwargs)  
    
    else:
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
        return {}  # Ensure it is always a dictionary
    
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
        unique_ids_str = kwargs.pop("unique_ids")
        # Skip if empty string
        if unique_ids_str:
            kwargs["unique_ids"] = unique_ids_str.split(",")

    if "unique_ids" in kwargs and len(kwargs) > 1:
        raise ValueError(
            "You can specify either unique_ids or metadata filters, not both"
        )

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
