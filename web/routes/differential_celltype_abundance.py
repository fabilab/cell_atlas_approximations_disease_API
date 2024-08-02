import os
import json
from flask import Blueprint, request, Response
from dotenv import load_dotenv
from google.cloud import storage
from utils.file_handling import run_function_flexible_infra
from utils.data_preprocessing import compute_diff_cell_abundance, get_metadata

# FIXME: what is this??
load_dotenv()

diff_celltype_abundance_bp = Blueprint('differential_celltype_abundance', __name__)
    

@diff_celltype_abundance_bp.route('/differential_cell_type_abundance', methods=['POST'])
def differential_cell_type_abundance():
    disease_keyword = request.args.get('disease_keyword', default='', type=str)
    id_list = request.args.get('unique_ids', default='', type=str) # comma-separated unique ids
    matching_datasets = get_metadata(disease_keyword, id_list.split(','))
    
    if len(matching_datasets) == 0:
        return []

    dataset_ids = [d['dataset_id'] for d in matching_datasets]
    diseases = [d['disease'] for d in matching_datasets]
    
    all_results = run_function_flexible_infra(compute_diff_cell_abundance)
    
    # Convert the list of OrderedDict to JSON string to preserve the order
    response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
    return Response(response_json, mimetype='application/json')
