"""
Create the obs database file with the metadata of all datasets.
"""
import os
import sys
import glob
import pathlib
import h5py
import hdf5plugin
import pandas as pd
import hashlib


def generate_unique_id(row):
    """Gnerate SHA256 hash for a given row"""
    row_string = ",".join(map(str, row))  # Convert row values to a single string
    m = hashlib.md5()
    m.update(row_string.encode("utf-8"))  # Update the hash with encoded row string
    return m.hexdigest() 


if __name__ == "__main__":

    data_folder = pathlib.Path(__file__).parent.absolute() / "static" / "atlas_data"

    h5_files = glob.glob(str(data_folder / "*.h5"))
    result = []
    for h5_fn in h5_files:
        dataset_id = os.path.basename(h5_fn).split(".")[0]
        with h5py.File(h5_fn, "r") as h5data:
            me = h5data["measurements"]
            for key in ["gene_expression", "chromatin_accessibility", "spatial"]:
                if key not in me:
                    continue
                group = me[key]
                col_names = [
                    group["groupby"]["names"].attrs[str(i)]
                    for i in range(group["groupby"].attrs["n_levels"])
                ]
                col_dtypes = [
                    group["groupby"]["dtypes"].attrs[str(i)]
                    for i in range(group["groupby"].attrs["n_levels"])
                ]
                res = {}
                for col_name, col_dtype in zip(col_names, col_dtypes):
                    # FIXME: they are not all strings
                    res[col_name] = group["obs"][col_name].asstr()[:]
                res["cell_count"] = group["obs"]["cell_count"][:]
                res = pd.DataFrame(res)
                res["dataset_id"] = dataset_id
                res = res[["dataset_id"] + col_names + ["cell_count"]]
                result.append(res)
    result = pd.concat(result)
    result.reset_index(drop=True, inplace=True)

    result["unique_id"] = result.apply(generate_unique_id, axis=1)

    result.to_csv(data_folder / "approximation_obs.csv", index=False)
