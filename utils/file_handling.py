import os
from dotenv import load_dotenv
from google.cloud import storage
import scquill

load_dotenv()

def download_blob(bucket_name, source_blob_name, destination_file_name):
    #  Downloads a blob from the bucket if it does not already exist locally.
    if not os.path.exists(destination_file_name):
        storage_client = storage.Client(project=os.getenv('GOOGLE_CLOUD_PROJECT'))
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(source_blob_name)
        blob.download_to_filename(destination_file_name)

def process_h5_file(file_path, keyword, compute_func, *args):
    # Assuming `file_path` is the path in the cloud storage bucket
    bucket_name = os.getenv('GOOGLE_CLOUD_BUCKET')
    local_path = f'/tmp/{os.path.basename(file_path)}'
    download_blob(bucket_name, file_path, local_path)
    
    # Continue processing the local file
    app = scquill.Approximation()
    app = app.read_h5(local_path)
    adata = app.to_anndata(groupby=('cell_type', 'disease'))
    
    if adata.obs[adata.obs['disease'].str.contains(keyword, case=False)].empty:
        return None
    
    dataset_id = os.path.basename(file_path).replace('.h5', '')
    return compute_func(adata, keyword, dataset_id, *args)