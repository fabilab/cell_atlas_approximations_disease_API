from flask import request, abort
from flask_restful import Resource

from api.v1.exceptions import model_exceptions
from api.v1.utils import (
    get_filter_kwargs,
    get_groupby_args,
)

from models.differential_gene_expression import get_diff_expression


class DifferentialGeneExpression(Resource):
    """Get differential gene expression base on user filter"""

    @model_exceptions
    def post(self):
        args = request.args

        differential_axis = args.get("differential_axis", "disease", type=str)
        groupby = get_groupby_args(args, ["tissue"])
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
        
        if "sex" in filters and "sex" not in groupby:
            groupby = ["sex"] + groupby
        
        feature = args.get("feature", None, type=str)
        number = int(request.args.get("top_n", default=0, type=int))
        if feature is not None and number > 0:
            abort(400, "Either feature or number must be provided, not both")
        
        if number == 0:
            number = 10

        method = args.get("method", "delta_fraction", type=str)

        diff_exp = get_diff_expression(
            differential_axis=differential_axis,
            groupby=groupby,
            number=number,
            feature=feature,
            method=method,
            **filters,
        )
        if len(diff_exp) == 0:
            return {
                "message": "No datasets found satisfying disease and cell type filter"
            }

        result = diff_exp.to_dict(orient="records")
        return result
