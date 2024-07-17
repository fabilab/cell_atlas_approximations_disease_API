from flask import Blueprint, request, Response
from dotenv import load_dotenv
from google.cloud import storage
from utils.file_handling import process_h5_file
from utils.data_preprocessing import compute_diff_expression
import logging
import json
import os

load_dotenv()

diff_gene_expression_bp = Blueprint('differential_gene_expression', __name__)

@diff_gene_expression_bp.route('/differential_gene_expression', methods=['POST'])
def differential_gene_expression():
    disease_keyword = request.args.get('disease_keyword', default='', type=str)
    id_list = request.args.get('unique_ids', default='', type=str) # comma-separated unique ids
    cell_type_keyword = request.args.get('cell_type', default='', type=str)
    top_n = int(request.args.get('top_n', default=10, type=int))
    h5_files_directory = 'compressed_data/h_sapiens/'  # Directory in the cloud storage bucket
    
    all_results = []
    # List all files in the bucket directory
    storage_client = storage.Client(project=os.getenv('GOOGLE_CLOUD_PROJECT'))
    blobs = storage_client.list_blobs(os.getenv('GOOGLE_CLOUD_BUCKET'), prefix=h5_files_directory)
    
    for blob in blobs:
        if blob.name.endswith('.h5'):
            result = process_h5_file(blob.name, disease_keyword, id_list, compute_diff_expression, top_n, cell_type_keyword)
            if result is not None:
                all_results.extend(result)

    # return jsonify(all_results)
    response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
    return Response(response_json, mimetype='application/json')