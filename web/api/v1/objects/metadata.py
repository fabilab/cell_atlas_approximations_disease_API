from flask import request
from flask_restful import Resource

from models.metadata import get_metadata

from api.v1.exceptions import model_exceptions
from api.v1.utils import get_filter_kwargs, validate_param_names


class Metadata(Resource):
    """
    Get a list of metadata based on disease, cell type, or tissue keywords
    This quickly help answer question such as:
    - is  there any record of the disease I am interested in?
    - What disease and datasets will contains the my cell type of interes?
    - What disease and datasets will contains the tissue I want?
    """

    @model_exceptions
    def get(self):
        args = request.args

        # Validate query parameter names to catch typos like "cell_type --> cell_typp"
        allowed_param_names = {
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
            ],
        )
        matching_metadata = get_metadata(**filters)

        # Convert the list of OrderedDict to JSON string to preserve the order
        result = matching_metadata.to_dict(orient="records")
        return result
