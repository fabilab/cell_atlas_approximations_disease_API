from flask import Blueprint, request, jsonify, Response
from utils.file_handling import process_h5_file
from utils.data_preprocessing import compute_diff_expression
import json
import os

# tutorial of Blueprint in flask: https://flask.palletsprojects.com/en/2.3.x/tutorial/views/
diff_gene_expression_bp = Blueprint('differential_gene_expression', __name__)

@diff_gene_expression_bp.route('/differential_gene_expression', methods=['GET'])
def differential_gene_expression():
    disease_keyword = request.args.get('keyword', default='', type=str)
    top_N = request.args.get('top_N', default=10, type=int)
    cell_type_keyword = request.args.get('celltype', default='', type=str)
    h5_files_directory = './data/h5_files_human'
    
    all_results = []
    for file_name in os.listdir(h5_files_directory):
        if file_name.endswith('.h5'):
            file_path = os.path.join(h5_files_directory, file_name)
            result = process_h5_file(
                file_path, 
                disease_keyword, 
                compute_diff_expression,
                top_N,
                cell_type_keyword
            )
            if result is not None:
                all_results.extend(result)
    

    response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
    return Response(response_json, mimetype='application/json')
