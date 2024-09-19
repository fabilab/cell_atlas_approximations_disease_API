"""Main module for API objects and endpoints"""
from api.v1 import api_dict as api_dict_v1


__all__ = (
    "api_dict",
)

api_dict = {
    "v1": api_dict_v1,
}