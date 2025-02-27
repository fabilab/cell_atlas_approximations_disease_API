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
    
    
    Optional Query Parameters:
    - `include_normal` (boolean, default=False): If `true`, includes the corresponding normal condition 
      when querying a disease. Only applicable when `disease` is provided.
    """

    @model_exceptions
    @required_parameters("features")
    def post(self):
        args = request.args
        
        # Optional argument to include the average expression of the normal condition.
        # This is only applicable when a disease condition is specified by the user.       
        include_normal = args.get("include_normal", "").lower() == "true"

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

        groupby = get_groupby_args(args, ["tissue", "disease"])
        feature_string = args.get("features", type=str)
        features = clean_feature_string(feature_string)
        
        result = get_average(
            features,
            groupby=groupby,
            include_normal=include_normal,
            **filters,
        )
        result = result.to_dict(orient="records")

        return result
