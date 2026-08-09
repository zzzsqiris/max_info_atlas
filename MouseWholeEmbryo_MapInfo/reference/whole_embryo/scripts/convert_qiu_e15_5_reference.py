#!/usr/bin/env python3
"""Convert the Qiu et al. run_20 E15.5 nuclei into a sparse AnnData file."""

import argparse
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


TASK_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = TASK_DIR / "data/raw/GSE228590"

CELL_ANNOTATION = RAW_DIR / "annotations/GSE228590_cell_annotation.run_20.csv.gz"
GENE_ANNOTATION = RAW_DIR / "annotations/GSE228590_gene_annotation.csv.gz"
COUNT_MATRIX = RAW_DIR / "counts/GSE228590_gene_count.run_20.mtx.gz"
OUTPUT = TASK_DIR / "data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad"
TEMP_OUTPUT = OUTPUT.with_suffix(".tmp.h5ad")

TARGET_DAY = "E15.5"
EXPECTED_SOURCE_GENES = 49_585
EXPECTED_SOURCE_CELLS = 827_534
EXPECTED_SOURCE_NNZ = 1_573_771_289
EXPECTED_TARGET_CELLS = 215_284
EXPECTED_EMBRYO = "embryo_59"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate metadata and the Matrix Market header without loading counts.",
    )
    return parser.parse_args()


def read_matrix_header():
    with gzip.open(COUNT_MATRIX, "rt") as handle:
        banner = handle.readline().strip()
        require(
            banner == "%%MatrixMarket matrix coordinate integer general",
            f"Unexpected Matrix Market banner: {banner}",
        )
        for line in handle:
            if not line.startswith("%"):
                shape = tuple(int(value) for value in line.split())
                break
        else:
            raise ValueError("Matrix Market size line was not found.")

    expected = (EXPECTED_SOURCE_GENES, EXPECTED_SOURCE_CELLS, EXPECTED_SOURCE_NNZ)
    require(shape == expected, f"Unexpected Matrix Market header: {shape}.")
    return shape


def read_metadata():
    cells = pd.read_csv(CELL_ANNOTATION, dtype=str)
    genes = pd.read_csv(GENE_ANNOTATION, dtype=str)

    require(len(cells) == EXPECTED_SOURCE_CELLS, "Unexpected run_20 cell count.")
    require(len(genes) == EXPECTED_SOURCE_GENES, "Unexpected gene count.")
    require(cells["cell_id"].is_unique, "Cell IDs are not unique.")
    require(genes["gene_ID"].is_unique, "Gene IDs are not unique.")
    require(not cells.isna().any().any(), "Cell annotation contains missing values.")
    require(
        not genes[["gene_ID", "gene_short_name"]].isna().any().any(),
        "Required gene annotation contains missing values.",
    )

    target_mask = cells["day"].eq(TARGET_DAY).to_numpy()
    require(target_mask.sum() == EXPECTED_TARGET_CELLS, "Unexpected E15.5 cell count.")
    target_embryos = cells.loc[target_mask, "embryo_id"].unique().tolist()
    require(target_embryos == [EXPECTED_EMBRYO], f"Unexpected E15.5 embryos: {target_embryos}")

    return cells, genes, target_mask


def collapse_gene_symbols(gene_by_cell, genes):
    """Sum source gene rows that have the same gene_short_name."""
    codes, gene_names = pd.factorize(genes["gene_short_name"], sort=False)
    require(np.all(codes >= 0), "A gene symbol could not be encoded.")

    source = gene_by_cell.tocoo(copy=False)
    collapsed = sparse.csc_matrix(
        (source.data, (codes[source.row], source.col)),
        shape=(len(gene_names), source.shape[1]),
        dtype=np.int32,
    )
    collapsed.sum_duplicates()

    grouped_ids = (
        genes.groupby("gene_short_name", sort=False, observed=True)["gene_ID"]
        .agg(";".join)
        .reindex(gene_names)
    )
    source_counts = (
        genes["gene_short_name"].value_counts(sort=False).reindex(gene_names)
    )
    var = pd.DataFrame(index=pd.Index(gene_names, name="gene_name"))
    var["gene_id"] = grouped_ids.to_numpy()
    var["source_feature_count"] = source_counts.to_numpy(dtype=np.int16)
    var["feature_type"] = "Gene Expression"
    return collapsed, var


def verify_output(path, expected_obs_names, expected_var_names):
    saved = ad.read_h5ad(path, backed="r")
    try:
        require(saved.shape == (EXPECTED_TARGET_CELLS, len(expected_var_names)), "Wrong output shape.")
        require(saved.obs_names.equals(expected_obs_names), "Saved cell order changed.")
        require(saved.var_names.equals(expected_var_names), "Saved gene order changed.")
        require(saved.obs_names.is_unique, "Saved cell IDs are not unique.")
        require(saved.var_names.is_unique, "Saved gene names are not unique.")
        require(getattr(saved.X, "format", None) in {"csr", "csc"}, "Saved X is not sparse.")
    finally:
        saved.file.close()


def main():
    args = parse_args()
    if OUTPUT.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists: {OUTPUT}. Use --overwrite to replace it.")

    cells, genes, target_mask = read_metadata()
    matrix_header = read_matrix_header()
    print(f"source cells: {len(cells)}")
    print(f"source genes: {len(genes)}")
    print(f"selected {TARGET_DAY} cells: {target_mask.sum()}")
    print(f"matrix header: {matrix_header}")
    if args.check_only:
        print("metadata and Matrix Market header checks passed")
        return
    print("reading full run_20 sparse Matrix Market file")

    with gzip.open(COUNT_MATRIX, "rb") as handle:
        gene_by_cell = sparse.coo_matrix(mmread(handle))

    require(
        gene_by_cell.shape == (EXPECTED_SOURCE_GENES, EXPECTED_SOURCE_CELLS),
        f"Unexpected matrix shape: {gene_by_cell.shape}.",
    )
    require(gene_by_cell.nnz == EXPECTED_SOURCE_NNZ, "Unexpected source nonzero count.")
    require(np.all(gene_by_cell.data >= 0), "Count matrix contains negative values.")
    require(
        gene_by_cell.data.max(initial=0) <= np.iinfo(np.int32).max,
        "A count exceeds int32.",
    )
    gene_by_cell.data = gene_by_cell.data.astype(np.int32, copy=False)

    print("converting source matrix to CSC and selecting E15.5 columns")
    gene_by_cell = gene_by_cell.tocsc()
    selected = gene_by_cell[:, target_mask]
    del gene_by_cell

    print("collapsing duplicated gene symbols")
    selected, var = collapse_gene_symbols(selected, genes)
    x = selected.T.tocsr(copy=False)
    x.sort_indices()
    del selected

    require(x.shape == (EXPECTED_TARGET_CELLS, len(var)), "Unexpected selected matrix shape.")
    require(np.all(x.data >= 0), "Selected matrix contains negative counts.")

    source_columns = np.flatnonzero(target_mask).astype(np.int32)
    obs = cells.loc[target_mask].copy()
    obs.index = pd.Index(obs.pop("cell_id"), name="cell_id")
    obs["source_matrix_column"] = source_columns
    obs["sample"] = pd.Categorical(["Qiu_E15.5_embryo_59"] * len(obs))
    obs["assay"] = pd.Categorical(["sci-RNA-seq3_nuclei"] * len(obs))

    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.uns["counts_state"] = "raw"
    adata.uns["source_accession"] = "GSE228590"
    adata.uns["source_run"] = "run_20"
    adata.uns["selected_day"] = TARGET_DAY
    adata.uns["source_matrix_shape"] = np.asarray(
        [EXPECTED_SOURCE_GENES, EXPECTED_SOURCE_CELLS], dtype=np.int64
    )
    adata.uns["duplicate_gene_symbols"] = (
        "Source feature rows with the same gene_short_name were summed."
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    TEMP_OUTPUT.unlink(missing_ok=True)
    print(f"writing temporary output: {TEMP_OUTPUT}")
    adata.write_h5ad(TEMP_OUTPUT)
    verify_output(TEMP_OUTPUT, adata.obs_names.copy(), adata.var_names.copy())
    TEMP_OUTPUT.replace(OUTPUT)
    print(f"final output: {OUTPUT}")
    print(f"verified shape: {adata.shape}")
    print(f"nonzero counts: {adata.X.nnz}")


if __name__ == "__main__":
    main()
