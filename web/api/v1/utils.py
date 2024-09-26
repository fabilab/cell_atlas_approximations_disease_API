"""Utility functions for the API."""


def get_optional_metadata_kwargs(args, columns):
    """Return a dictionary of obs metadata keyword arguments to be used as database filters."""
    kwargs = {}
    for column in columns:
        value = args.get(column, default="", type=str)
        if value != "":
            kwargs[column] = value
    return kwargs
