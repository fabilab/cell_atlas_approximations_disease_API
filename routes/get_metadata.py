from flask import Blueprint, request, Response
from dotenv import load_dotenv
from google.cloud import storage
from utils.file_handling import process_h5_file
from utils.data_preprocessing import get_metadata
import json
import os

load_dotenv()

get_metadata_bp = Blueprint('metadata', __name__)

@get_metadata_bp.route('/metadata', methods=['GET'])
def get_metadata_route():
    disease_keyword = request.args.get('disease_keyword', default='', type=str)
    h5_files_directory = 'compressed_data/h_sapiens/'  # Directory in the cloud storage bucket
    
    all_results = []
    
    # List all files in the bucket directory
    storage_client = storage.Client(project=os.getenv('GOOGLE_CLOUD_PROJECT'))
    blobs = storage_client.list_blobs(os.getenv('GOOGLE_CLOUD_BUCKET'), prefix=h5_files_directory)
    
    for blob in blobs:
        if blob.name.endswith('.h5'):
            result = process_h5_file(blob.name, disease_keyword, '', get_metadata)
            if result is not None:
                all_results.extend(result)
    
    # Convert the list of OrderedDict to JSON string to preserve the order
    response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
    return Response(response_json, mimetype='application/json')
