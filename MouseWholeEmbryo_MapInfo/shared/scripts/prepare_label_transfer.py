#!/usr/bin/env python3
"""Prepare shared-gene features and correlation-kNN chunks."""

import argparse
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse

from max_info_atlases.uge.job_generator import create_chunk_files, write_job_list


PIPELINE_DIR = Path(__file__).resolve().parents[2]
DEFAULT_REFERENCE = (
    PIPELINE_DIR
    / "reference/whole_embryo/data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad"
)
DEFAULT_TARGET = (
    PIPELINE_DIR
    / "datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad"
)
DEFAULT_OUTPUT_DIR = PIPELINE_DIR / "datasets/visium_hd/results/label_transfer/work"
EXPECTED_REFERENCE_CELLS = 215_284
EXPECTED_TARGET_CELLS = 449_700
EXPECTED_SHARED_GENES = 0


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--target-cells-per-job", type=int, default=500)
    parser.add_argument(
        "--expected-reference-cells",
        type=int,
        default=EXPECTED_REFERENCE_CELLS,
    )
    parser.add_argument(
        "--expected-target-cells",
        type=int,
        default=EXPECTED_TARGET_CELLS,
    )
    parser.add_argument(
        "--expected-shared-genes",
        type=int,
        default=EXPECTED_SHARED_GENES,
        help="Set to 0 to record the observed intersection without a fixed count.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def to_dense_float64(matrix):
    if sparse.issparse(matrix):
        return matrix.toarray().astype(np.float64, copy=False)
    return np.asarray(matrix, dtype=np.float64)


def save_array(path, values):
    np.save(path, values, allow_pickle=False)


def main():
    args = parse_args()
    if args.target_cells_per_job < 1:
        raise ValueError("--target-cells-per-job must be positive.")
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(
            f"Output directory is not empty: {args.output_dir}. Use --overwrite."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    chunks_dir = args.output_dir / "chunks"
    if chunks_dir.exists() and args.overwrite:
        shutil.rmtree(chunks_dir)

    reference = ad.read_h5ad(args.reference, backed="r")
    target = ad.read_h5ad(args.target, backed="r")
    try:
        if reference.n_obs != args.expected_reference_cells:
            raise ValueError(
                f"Expected {args.expected_reference_cells} reference cells, "
                f"found {reference.n_obs}."
            )
        if target.n_obs != args.expected_target_cells:
            raise ValueError(
                f"Expected {args.expected_target_cells} target cells, found {target.n_obs}."
            )
        if not reference.obs_names.is_unique or not target.obs_names.is_unique:
            raise ValueError("Reference and target cell IDs must be unique.")
        if not reference.var_names.is_unique or not target.var_names.is_unique:
            raise ValueError("Reference and target gene names must be unique.")

        reference_genes = set(reference.var_names)
        shared_genes = [gene for gene in target.var_names if gene in reference_genes]
        if args.expected_shared_genes > 0 and len(shared_genes) != args.expected_shared_genes:
            raise ValueError(
                f"Expected {args.expected_shared_genes} shared genes, found {len(shared_genes)}."
            )

        reference_shared = reference[:, shared_genes].to_memory()
        target_shared = target[:, shared_genes].to_memory()
        reference_names = reference.obs_names.to_numpy(dtype=str)
        target_names = target.obs_names.to_numpy(dtype=str)
    finally:
        reference.file.close()
        target.file.close()

    if not reference_shared.var_names.equals(target_shared.var_names):
        raise ValueError("Reference and target shared-gene orders differ.")

    # Match the established LabelTransfer pipeline exactly.
    sc.pp.normalize_total(reference_shared, target_sum=10_000)
    sc.pp.log1p(reference_shared)
    sc.pp.normalize_total(target_shared, target_sum=10_000)
    sc.pp.log1p(target_shared)

    reference_values = to_dense_float64(reference_shared.X)
    target_values = to_dense_float64(target_shared.X)
    if not np.isfinite(reference_values).all() or not np.isfinite(target_values).all():
        raise ValueError("Normalized shared-gene features contain non-finite values.")

    retained = (np.var(reference_values, axis=0) > 0) & (
        np.var(target_values, axis=0) > 0
    )
    retained_genes = np.asarray(shared_genes, dtype=str)[retained]
    reference_values = reference_values[:, retained]
    target_values = target_values[:, retained]
    if len(retained_genes) == 0:
        raise ValueError("No variable shared genes remain.")
    reference_cell_variance = np.var(reference_values, axis=1)
    reference_source_indices = np.flatnonzero(reference_cell_variance > 0)
    excluded_reference_cells = len(reference_values) - len(reference_source_indices)
    reference_values = reference_values[reference_source_indices]
    reference_names = reference_names[reference_source_indices]
    if len(reference_values) == 0:
        raise ValueError("No usable reference cells remain for correlation kNN.")
    target_variance = np.var(target_values, axis=1)
    target_source_indices = np.flatnonzero(target_variance > 0)
    excluded_zero_variance_target_cells = len(target_values) - len(target_source_indices)
    if len(target_source_indices) == 0:
        raise ValueError("No target cells have usable variance across retained genes.")
    target_values = target_values[target_source_indices]
    target_names = target_names[target_source_indices]

    save_array(args.output_dir / "reference_features.npy", reference_values)
    save_array(args.output_dir / "target_features.npy", target_values)
    save_array(args.output_dir / "reference_obs_names.npy", reference_names)
    save_array(
        args.output_dir / "reference_source_indices.npy",
        reference_source_indices.astype(np.int32, copy=False),
    )
    save_array(
        args.output_dir / "target_source_indices.npy",
        target_source_indices.astype(np.int32, copy=False),
    )
    save_array(args.output_dir / "target_obs_names.npy", target_names)
    save_array(args.output_dir / "retained_genes.npy", retained_genes)

    jobs = []
    for start in range(0, len(target_names), args.target_cells_per_job):
        stop = min(start + args.target_cells_per_job, len(target_names))
        jobs.append(f"{start}\t{stop}")

    jobs_file = args.output_dir / "jobs.txt"
    write_job_list(jobs, jobs_file)
    chunk_files, _ = create_chunk_files(
        job_list=jobs,
        output_dir=chunks_dir,
        chunk_size=1,
    )

    manifest = {
        "reference": str(args.reference),
        "target": str(args.target),
        "source_reference_cells": len(reference_source_indices) + excluded_reference_cells,
        "usable_reference_cells": len(reference_names),
        "excluded_zero_variance_reference_cells": excluded_reference_cells,
        "target_cells": args.expected_target_cells,
        "target_cells_for_knn": len(target_names),
        "excluded_zero_variance_target_cells": excluded_zero_variance_target_cells,
        "initial_shared_genes": len(shared_genes),
        "retained_shared_genes": len(retained_genes),
        "target_cells_per_job": args.target_cells_per_job,
        "jobs": len(jobs),
        "normalization": "normalize_total(target_sum=10000), log1p",
        "distance": "correlation",
    }
    with open(args.output_dir / "manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    print(f"reference_shape: {reference_values.shape}")
    print(f"excluded_zero_variance_reference_cells: {excluded_reference_cells}")
    print(f"target_shape: {target_values.shape}")
    print(f"excluded_zero_variance_target_cells: {excluded_zero_variance_target_cells}")
    print(f"retained_shared_genes: {len(retained_genes)}")
    print(f"jobs: {len(jobs)}")
    print(f"chunk_files: {len(chunk_files)}")
    print(f"output_dir: {args.output_dir}")


if __name__ == "__main__":
    main()
