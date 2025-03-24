from flask import request
from flask_restful import Resource

from api.v1.exceptions import model_exceptions
from api.v1.utils import (
    get_filter_kwargs,
    get_groupby_args,
    validate_param_names,
)

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

        differential_axis = args.get("differential_axis", "disease", type=str)
        groupby = get_groupby_args(args, ["tissue_general"])
        
        # Validate query parameter names to catch typos like "cell_type --> cell_typp"
        allowed_param_names = {
            "differential_axis",
            "disease",
            "cell_type",
            "tissue",
            "sex",
            "development_stage",
        }

        validate_param_names(args, allowed_param_names)

        filters = get_filter_kwargs(
            args,
            [
                "disease",
                "cell_type",
                "tissue",
                "sex",
                "development_stage",
                # "unique_ids",
            ],
        )

        diff_cell_abundance = get_diff_cell_abundance(
            differential_axis=differential_axis,
            groupby=groupby,
            **filters,
        )
        if len(diff_cell_abundance) == 0:
            filt_string = ", ".join(filters.keys())
            return {
                "message": f"No datasets found satisfying the requested filters: {filt_string}"
            }

        result = diff_cell_abundance.to_dict(orient="records")
        return result
