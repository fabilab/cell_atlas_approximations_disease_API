from flask import request
from flask_restful import Resource

from models import (
    get_highest_measurement,
)

from api.v1.exceptions import (
    model_exceptions,
    required_parameters,
)
from api.v1.utils import get_groupby_args, validate_param_names


class HighestMeasurement(Resource):
    """
    API Resource class to handle requests for the highest expression measurement
    of a given feature (gene) across multiple all diseases and their datasets . It returns the top n
    disease and cell types combination with the highest expression of the specified gene.
    """

    @model_exceptions
    @required_parameters("feature")
    def post(self):
        args = request.args

        groupby = get_groupby_args(args, ["tissue", "disease"])
        feature = args.get("feature")
        number = args.get("number", default=10)
        
        # Validate query parameter names to catch typos like "feature --> feture"
        allowed_param_names = {
            "feature",
            "number",
        }

        validate_param_names(args, allowed_param_names)

        # Get the top N highest expressors from all the expression results
        highest_measurement = get_highest_measurement(
            feature,
            number,
            groupby=groupby,
        )
        result = highest_measurement.to_dict(orient="records")

        return result
