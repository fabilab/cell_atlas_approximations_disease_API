import json
import os
from config import configuration as config
from flask import Response, request
from flask_restful import Resource, abort

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

        number = int(request.args.get("top_n", default=10, type=int))

        matching_datasets = get_metadata(
            filters["disease"], filters["cell_type"], filters["unique_ids"]
        )

        diff_exp = get_diff_expression(number=number, **filters)
        if len(diff_exp) == 0:
            return {
                "message": "No datasets found satisfying disease and cell type filter"
            }

        # return jsonify(all_results)
        response_json = json.dumps(diff_exp, ensure_ascii=False, indent=4)
        return Response(response_json, mimetype="application/json")
