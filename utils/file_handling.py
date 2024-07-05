import scquill

# Initialize the scquill Approximation object and read the HDF5 file to obtain the AnnData object
# https://github.com/fabilab/scquill/tree/main
def process_h5_file(file_path, keyword, compute_func, *args):
    print(file_path)
    app = scquill.Approximation()
    app = app.read_h5(file_path)
    adata = app.to_anndata(groupby=('cell_type', 'disease'))
    
    # filter out datasets that don't contain the user keyword in the disease column
    if adata.obs[adata.obs['disease'].str.contains(keyword, case=False)].empty:
        return None
    
    dataset_id = file_path.split('/')[-1].replace('.h5', '')
    return compute_func(adata, keyword, dataset_id, *args)
