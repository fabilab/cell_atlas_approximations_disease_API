from flask import request
from flask_restful import Resource

from models import (
    get_highest_measurement,
)

from api.v1.exceptions import (
    model_exceptions,
    required_parameters,
)


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

        feature = args.get("feature")
        number = args.get("number", default=10)

        # Get the top N highest expressors from all the expression results
        highest_measurement = get_highest_measurement(feature, number)
        result = highest_measurement.to_dict(orient="records")

        return result
