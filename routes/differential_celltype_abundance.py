# routes/differential_celltype_abundance.py

from flask import Blueprint, request, jsonify, Response
import os
from utils.file_handling import process_h5_file
from utils.data_preprocessing import compute_diff_cell_abundance
import json

diff_celltype_abundance_bp = Blueprint('differential_celltype_abundance', __name__)

@diff_celltype_abundance_bp.route('/differential_cell_type_abundance', methods=['GET'])
def differential_cell_type_abundance():
    disease_keyword = request.args.get('keyword', default='', type=str)
    h5_files_directory = './data/h5_files_human'
    
    all_results = []
    for file_name in os.listdir(h5_files_directory):
        if file_name.endswith('.h5'):
            file_path = os.path.join(h5_files_directory, file_name)
            result = process_h5_file(file_path, disease_keyword, compute_diff_cell_abundance)
            if result is not None:
                all_results.extend(result)
    
    # return jsonify(all_results)
    # Convert the list of OrderedDict to JSON string to preserve the order
    response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
    return Response(response_json, mimetype='application/json')
