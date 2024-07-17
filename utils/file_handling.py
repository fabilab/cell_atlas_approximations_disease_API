import os
from dotenv import load_dotenv
from google.cloud import storage
import scquill
import h5py
from inspect import signature

load_dotenv()

def download_blob(bucket_name, source_blob_name, destination_file_name):
    #  Downloads a blob from the bucket if it does not already exist locally.
    if not os.path.exists(destination_file_name):
        storage_client = storage.Client(project=os.getenv('GOOGLE_CLOUD_PROJECT'))
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(source_blob_name)
        blob.download_to_filename(destination_file_name)

def process_h5_file(file_path, disease_keyword, unique_id_list, compute_func, *args):
    bucket_name = os.getenv('GOOGLE_CLOUD_BUCKET')
    local_path = f'/tmp/{os.path.basename(file_path)}'
    download_blob(bucket_name, file_path, local_path)
    
    app = scquill.Approximation()
    app = app.read_h5(local_path)
    adata = app.to_anndata(groupby=('cell_type', 'disease'))
    
    with h5py.File(local_path, 'r') as compressed_h5:
            
        metadata = {
            'unique_ids': compressed_h5.attrs['ids'],
            'diseases': compressed_h5.attrs['diseases'],
            'dataset_id': compressed_h5.attrs['dataset_id'],
            'cell_types': compressed_h5.attrs['cell_types'],
            'unit': compressed_h5.attrs['unit'],
            'log_transformed': compressed_h5.attrs['log_transformed'],
            'has_normal_baseline': compressed_h5.attrs['has_normal_baseline'],
        }

        # Create a dictionary pairing unique IDs with diseases
        id_disease_pair = {unique_id.decode('utf-8'): disease.decode('utf-8') for unique_id, disease in zip(metadata['unique_ids'], metadata['diseases'])}
    
    unique_ids = [item for item in unique_id_list.split(',') if item != '']

    if adata.obs[adata.obs['disease'].str.contains(disease_keyword, case=False)].empty:
        return None

    # if user supplied a keyword to filter on disease
    # check if any of the disease in metadata contains the keyword as a substring
    # if None, then return
    # if len([encoded_d for encoded_d in metadata['diseases'] if disease_keyword in encoded_d.decode()]) == 0:
    #     return None
    
    # if user supplied unique ids
    if len(unique_ids) > 0:
        # check if any of the user supplied unique ids can be found in this dataset
        unique_ids_found_in_dataset = [unique_id for unique_id in unique_ids if unique_id in [uid.decode() for uid in metadata['unique_ids']]]
        if len(unique_ids_found_in_dataset) == 0:
            return None
        # here, we need to extract the disease keyword using the unique id and pass it to data preprocessing function
        if not disease_keyword:
            disease_keyword = id_disease_pair[unique_ids_found_in_dataset[0]]
        
    dataset_id = os.path.basename(file_path).replace('.h5', '')
    
    # TO FIX: hard coded here, might not work in the future. 
    # Check the number of arguments that the compute function expects
    sig = signature(compute_func)
    param_count = len(sig.parameters)
    
    if param_count == 2:
        return compute_func(disease_keyword, metadata)
    elif param_count == 3:
        return compute_func(adata, disease_keyword, dataset_id)
    else:
        return compute_func(adata, disease_keyword, dataset_id, metadata, *args)