from flask import request, abort
from flask_restful import Resource

from api.v1.exceptions import model_exceptions
from api.v1.utils import get_optional_metadata_kwargs

from models.differential_gene_expression import get_diff_expression


class DifferentialGeneExpression(Resource):
    """Get differential gene expression base on user filter"""

    @model_exceptions
    def post(self):
        args = request.args

        filters = get_optional_metadata_kwargs(
            args, ["disease", "cell_type", "tissue", "sex", "development_stage"]
        )
        unique_ids_str = request.args.get("unique_ids", default="", type=str)
        if unique_ids_str:
            filters["unique_ids"] = unique_ids_str.split(",")

        feature = args.get("feature", "", type=str)
        number = int(request.args.get("top_n", default=0, type=int))
        if feature is not None and number > 0:
            abort(400, "Either feature or number must be provided, not both")
        elif feature is None:
            number = 10

        diff_exp = get_diff_expression(number=number, feature=feature, **filters)
        if len(diff_exp) == 0:
            return {
                "message": "No datasets found satisfying disease and cell type filter"
            }

        result = diff_exp.to_dict(orient="records")
        return result
