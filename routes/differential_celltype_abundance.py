from flask import Blueprint, request, Response
from dotenv import load_dotenv
from google.cloud import storage
from utils.file_handling import process_h5_file
from utils.data_preprocessing import compute_diff_cell_abundance
import json
import os

load_dotenv()

diff_celltype_abundance_bp = Blueprint('differential_celltype_abundance', __name__)

@diff_celltype_abundance_bp.route('/differential_cell_type_abundance', methods=['POST'])
def differential_cell_type_abundance():
    disease_keyword = request.args.get('disease_keyword', default='', type=str)
    id_list = request.args.get('unique_ids', default='', type=str) # comma-separated unique ids
    h5_files_directory = 'compressed_data/h_sapiens/'  # Directory in the cloud storage bucket
    
    all_results = []
    
    # List all files in the bucket directory
    storage_client = storage.Client(project=os.getenv('GOOGLE_CLOUD_PROJECT'))
    blobs = storage_client.list_blobs(os.getenv('GOOGLE_CLOUD_BUCKET'), prefix=h5_files_directory)
    
    for blob in blobs:
        if blob.name.endswith('.h5'):
            result = process_h5_file(blob.name, disease_keyword, id_list, compute_diff_cell_abundance)
            if result is not None:
                all_results.extend(result)
    
    
    # return jsonify(all_results)
    # Convert the list of OrderedDict to JSON string to preserve the order
    response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
    return Response(response_json, mimetype='application/json')