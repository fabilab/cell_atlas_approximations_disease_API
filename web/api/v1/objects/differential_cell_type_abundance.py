from flask_restful import Resource

from api.v1.exceptions import model_exceptions
from api.v1.utils import get_optional_metadata_kwargs

from models import (
    get_diff_cell_abundance,
)


class DifferentialCellTypeAbundance(Resource):
    """
    Get differential cell type abundance for a disease or tissue of interest across multiple datasets.

    This API helps answer questions such as:
    1. "In COVID-19, what is the differential cell abundance of each cell type between normal and disease states?"
    2. "Across all diseases affecting the lung, what is the differential cell abundance of each cell type in that tissue between normal and disease states?"
    """

    @model_exceptions
    def post(self):
        args = request.args

        filters = get_optional_metadata_kwargs(
            args, ["disease", "cell_type", "tissue", "sex"]
        )
        unique_ids_str = request.args.get("unique_ids", default="", type=str)
        if unique_ids_str:
            filters["unique_ids"] = unique_ids_str.split(",")

        diff_cell_abundance = get_diff_cell_abundance(**filters)
        if len(diff_cell) == 0:
            filt_string = ", ".join(filters.keys())
            return {
                "message": f"No datasets found satisfying the requested filters: {filt_string}"
            }

        result = diff_cell_abundance.to_dict(orient="records")
        return result
