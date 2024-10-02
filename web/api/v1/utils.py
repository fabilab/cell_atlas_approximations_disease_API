"""Utility functions for the API."""


# FIXME: there is a logical fallacy in allowing both unique_ids and metadta filters, deal with it at some point
def get_filter_kwargs(args, columns):
    """Return a dictionary of obs metadata keyword arguments to be used as database filters."""
    kwargs = {}
    for column in columns:
        value = args.get(column, default="", type=str)
        if value != "":
            kwargs[column] = value

    kwargs = _clean_metadata_kwargs(kwargs)

    return kwargs


def _clean_metadata_kwargs(kwargs):
    """Clean metadata filters"""
    if "sex" in kwargs:
        value = _normalise_sex_string(kwargs.pop("sex"))
        if value is not None:
            kwargs["sex"] = value

    # NOTE: This is a temporary fix re renaming of column
    if "development_stage" in kwargs:
        kwargs["development_stage_general"] = kwargs.pop("development_stage")

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
