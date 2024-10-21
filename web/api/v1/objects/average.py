from flask import request
from flask_restful import Resource

from models import (
    get_average,
)

from api.v1.utils import (
    clean_feature_string,
    get_filter_kwargs,
)
from api.v1.exceptions import (
    model_exceptions,
    required_parameters,
)
from api.v1.utils import get_groupby_args


class Average(Resource):
    """
    API Resource class to handle requests for the average measurement
    of select features (gene) across multiple all diseases and their datasets.
    """

    @model_exceptions
    @required_parameters("features")
    def post(self):
        args = request.args

        groupby = get_groupby_args(args, ["tissue", "disease"])
        feature_string = args.get("features", type=str)
        features = clean_feature_string(feature_string)

        filters = get_filter_kwargs(
            args,
            [
                "disease",
                "cell_type",
                "tissue",
                "sex",
                "development_stage",
                "unique_ids",
            ],
        )

        result = get_average(
            features,
            groupby=groupby,
            **filters,
        )
        result = result.to_dict(orient="records")

        return result
