#!/usr/bin/env python3
"""Convert the 10x Visium HD segmented output to the Map Information AnnData format."""

import argparse
import shutil
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse


TASK_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = TASK_DIR / "data/raw"
OUTPUT = TASK_DIR / "data/processed/mouse_embryo_visium_hd_raw.h5ad"
MATRIX_TAR = RAW_DIR / "Visium_HD_3prime_Mouse_Embryo_segmented_outputs.tar.gz"
MAPPING = RAW_DIR / "Visium_HD_3prime_Mouse_Embryo_barcode_mappings.parquet"
MATRIX_MEMBER = "segmented_outputs/raw_feature_cell_matrix.h5"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def decode(values):
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values],
        dtype=str,
    )


def extract_matrix_h5(destination):
    with tarfile.open(MATRIX_TAR, mode="r:gz") as archive:
        member = archive.getmember(MATRIX_MEMBER)
        source = archive.extractfile(member)
        require(source is not None, f"Could not read {MATRIX_MEMBER} from {MATRIX_TAR}.")
        with source, destination.open("wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)


def read_matrix(path):
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(value) for value in matrix["shape"][:])
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        barcodes = pd.Index(decode(matrix["barcodes"][:]), name="cell_id")
        features = pd.DataFrame(
            {
                "gene_id": decode(matrix["features"]["id"][:]),
                "gene_name": decode(matrix["features"]["name"][:]),
                "feature_type": decode(matrix["features"]["feature_type"][:]),
            }
        )

    require(shape == (len(features), len(barcodes)), f"Matrix shape does not match metadata: {shape}.")
    require(barcodes.is_unique, "Visium cell barcodes are not unique.")
    require(features["feature_type"].eq("Gene Expression").all(), "Non-gene features are present.")
    require(features["gene_name"].notna().all(), "Gene names contain missing values.")

    gene_by_cell = sparse.csc_matrix((data, indices, indptr), shape=shape)
    codes, gene_names = pd.factorize(features["gene_name"], sort=False)
    coo = gene_by_cell.tocoo(copy=False)
    collapsed = sparse.csr_matrix(
        (coo.data, (codes[coo.row], coo.col)),
        shape=(len(gene_names), gene_by_cell.shape[1]),
        dtype=np.int32,
    )
    collapsed.sum_duplicates()
    x = collapsed.T.tocsr(copy=False)
    x.sort_indices()

    grouped_ids = (
        features.groupby("gene_name", sort=False, observed=True)["gene_id"]
        .agg(";".join)
        .reindex(gene_names)
    )
    source_counts = features["gene_name"].value_counts(sort=False).reindex(gene_names)
    var = pd.DataFrame(index=pd.Index(gene_names, name="gene_name"))
    var["gene_id"] = grouped_ids.to_numpy()
    var["source_feature_count"] = source_counts.to_numpy(dtype=np.int16)
    var["feature_type"] = "Gene Expression"
    return x, barcodes, var


def read_coordinates(barcodes):
    mapping = pq.read_table(MAPPING, columns=["square_002um", "cell_id"]).to_pandas()
    mapping = mapping.loc[mapping["cell_id"].notna()].copy()
    parts = mapping["square_002um"].str.extract(r"^s_002um_(\d+)_(\d+)-1$")
    require(not parts.isna().any().any(), "Could not parse a 2um square barcode.")
    mapping["row"] = parts[0].astype(np.float64)
    mapping["col"] = parts[1].astype(np.float64)
    centroids = mapping.groupby("cell_id", sort=False)[["col", "row"]].mean()
    centroids = centroids.reindex(barcodes)
    require(not centroids.isna().any().any(), "Some matrix barcodes lack spatial coordinates.")
    xy = centroids.to_numpy(dtype=np.float64)
    require(np.isfinite(xy).all(), "Spatial coordinates contain non-finite values.")
    require(len(np.unique(xy, axis=0)) == len(xy), "Spatial coordinates contain duplicates.")
    return xy


def verify_output(path, expected_shape, expected_obs, expected_var):
    saved = ad.read_h5ad(path, backed="r")
    try:
        require(saved.shape == expected_shape, f"Wrong saved shape: {saved.shape}.")
        require(saved.obs_names.equals(expected_obs), "Saved cell order changed.")
        require(saved.var_names.equals(expected_var), "Saved gene order changed.")
        require(saved.obs_names.is_unique and saved.var_names.is_unique, "Saved IDs are not unique.")
        require(getattr(saved.X, "format", None) in {"csr", "csc"}, "Saved X is not sparse.")
        require("sample" in saved.obs, "Saved AnnData lacks obs['sample'].")
        require("spatial" in saved.obsm, "Saved AnnData lacks obsm['spatial'].")
        require(saved.obsm["spatial"].shape == (saved.n_obs, 2), "Wrong spatial coordinate shape.")
    finally:
        saved.file.close()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    require(MATRIX_TAR.exists(), f"Missing raw matrix archive: {MATRIX_TAR}")
    require(MAPPING.exists(), f"Missing barcode mapping: {MAPPING}")
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists: {args.output}. Use --overwrite.")

    with tempfile.TemporaryDirectory(prefix="visium_hd_") as temp_dir:
        matrix_h5 = Path(temp_dir) / "raw_feature_cell_matrix.h5"
        print(f"extracting matrix HDF5 to temporary storage: {matrix_h5}")
        extract_matrix_h5(matrix_h5)
        x, barcodes, var = read_matrix(matrix_h5)

    xy = read_coordinates(barcodes)
    obs = pd.DataFrame(index=barcodes)
    obs["sample"] = pd.Categorical(["MouseEmbryo_VisiumHD"] * len(barcodes))
    obs["assay"] = pd.Categorical(["Visium HD 3-prime Gene Expression"] * len(barcodes))

    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["spatial"] = xy
    adata.uns["counts_state"] = "raw"
    adata.uns["source_matrix"] = "10x segmented_outputs/raw_feature_cell_matrix.h5"
    adata.uns["spatial_coordinate_definition"] = (
        "Mean 2um square column,row for each segmented cell from barcode_mappings.parquet."
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = args.output.with_suffix(".tmp.h5ad")
    temporary_output.unlink(missing_ok=True)
    adata.write_h5ad(temporary_output)
    verify_output(temporary_output, adata.shape, adata.obs_names.copy(), adata.var_names.copy())
    temporary_output.replace(args.output)
    print(f"verified output: {args.output}")
    print(f"shape: {adata.shape}")
    print(f"nonzero counts: {adata.X.nnz}")


if __name__ == "__main__":
    main()
