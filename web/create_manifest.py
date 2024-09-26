import os
import sys
import glob
import pathlib
import h5py
import hdf5plugin
import pandas as pd


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
                res["ncells"] = group["obs"]["cell_count"][:]
                res = pd.DataFrame(res)[col_names]
                res["dataset_id"] = dataset_id
                result.append(res)
    result = pd.concat(result)

    result.to_csv(data_folder / "approximation_obs.csv", index=False)
