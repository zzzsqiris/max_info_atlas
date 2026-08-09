#!/usr/bin/env python3
"""Combine kNN chunks and transfer all reference labels to a spatial target."""

import argparse
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PIPELINE_DIR = Path(__file__).resolve().parents[2]
DEFAULT_TARGET = (
    PIPELINE_DIR
    / "datasets/visium_hd/data/processed/mouse_embryo_visium_hd_raw.h5ad"
)
DEFAULT_LABELS_DIR = (
    PIPELINE_DIR
    / "reference/whole_embryo/results/reference_clustering/clustering/LeidenPCA50Correlation"
)
VISIUM_DIR = PIPELINE_DIR / "datasets/visium_hd"
DEFAULT_WORK_DIR = VISIUM_DIR / "results/label_transfer/work"
DEFAULT_CHUNKS_DIR = VISIUM_DIR / "results/label_transfer/chunk_results"
DEFAULT_OUTPUT = VISIUM_DIR / "data/processed/mouse_embryo_visium_hd_label_transfer.h5ad"
DEFAULT_SUMMARY = VISIUM_DIR / "results/label_transfer/label_transfer_summary.csv"
EXPECTED_RESOLUTIONS = 50
EXPECTED_SOURCE_REFERENCE_CELLS = 215_284
EXPECTED_TARGET_CELLS = 449_700
REFERENCE_SECTION = "Qiu_E15.5_embryo_59"
UNASSIGNED_LABEL = "unassigned_zero_shared_expression"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--labels-dir", type=Path, default=DEFAULT_LABELS_DIR)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--chunks-dir", type=Path, default=DEFAULT_CHUNKS_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument(
        "--expected-source-reference-cells",
        type=int,
        default=EXPECTED_SOURCE_REFERENCE_CELLS,
    )
    parser.add_argument(
        "--expected-target-cells",
        type=int,
        default=EXPECTED_TARGET_CELLS,
    )
    parser.add_argument("--reference-section", default=REFERENCE_SECTION)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def resolution_value(path):
    match = re.fullmatch(r"res_(\d+)p(\d+)", path.parent.name)
    if not match:
        raise ValueError(f"Unexpected resolution directory: {path.parent.name}")
    return float(f"{match.group(1)}.{match.group(2)}")


def select_vote(neighbor_labels, neighbor_distances, weights):
    labels = np.unique(neighbor_labels)
    scores = np.asarray(
        [weights[neighbor_labels == label].sum() for label in labels],
        dtype=np.float64,
    )
    mean_distances = np.asarray(
        [neighbor_distances[neighbor_labels == label].mean() for label in labels],
        dtype=np.float64,
    )
    candidates = np.flatnonzero(scores == scores.max())
    if len(candidates) > 1:
        candidate_means = mean_distances[candidates]
        candidates = candidates[candidate_means == candidate_means.min()]
    winner = int(candidates[0])
    sorted_scores = np.sort(scores)
    second_score = sorted_scores[-2] if len(sorted_scores) > 1 else 0.0
    return labels[winner], scores[winner], second_score


def load_neighbors(chunks_dir, n_target, k):
    indices = np.empty((n_target, k), dtype=np.int32)
    distances = np.empty((n_target, k), dtype=np.float32)
    coverage = np.zeros(n_target, dtype=np.int8)

    files = sorted(chunks_dir.glob("neighbors_*.npz"))
    if not files:
        raise FileNotFoundError(f"No neighbor chunks found in {chunks_dir}.")

    for path in files:
        with np.load(path, allow_pickle=False) as chunk:
            start = int(chunk["start"])
            stop = int(chunk["stop"])
            chunk_indices = chunk["indices"]
            chunk_distances = chunk["distances"]
        expected_shape = (stop - start, k)
        if chunk_indices.shape != expected_shape or chunk_distances.shape != expected_shape:
            raise ValueError(f"Wrong neighbor shape in {path}.")
        if start < 0 or stop > n_target or stop <= start:
            raise ValueError(f"Invalid target range in {path}.")
        if np.any(coverage[start:stop]):
            raise ValueError(f"Overlapping target range in {path}.")
        indices[start:stop] = chunk_indices
        distances[start:stop] = chunk_distances
        coverage[start:stop] += 1

    if not np.all(coverage == 1):
        missing = int(np.count_nonzero(coverage == 0))
        raise ValueError(f"Neighbor chunks do not cover {missing} target cells.")
    return indices, distances, len(files)


def transfer_one_resolution(reference_labels, indices, distances):
    n_target, k = indices.shape
    unweighted = np.empty(n_target, dtype=reference_labels.dtype)
    weighted = np.empty(n_target, dtype=reference_labels.dtype)
    vote_fraction = np.empty(n_target, dtype=np.float32)
    weighted_proportion = np.empty(n_target, dtype=np.float32)
    weighted_margin = np.empty(n_target, dtype=np.float32)

    for cell_index in range(n_target):
        cell_distances = distances[cell_index].astype(np.float64, copy=False)
        cell_labels = reference_labels[indices[cell_index]]
        weights = 1.0 / (cell_distances + 1e-8)

        weighted_label, weighted_score, second_score = select_vote(
            cell_labels, cell_distances, weights
        )
        unweighted_label, unweighted_score, _ = select_vote(
            cell_labels,
            cell_distances,
            np.ones(k, dtype=np.float64),
        )
        weighted[cell_index] = weighted_label
        unweighted[cell_index] = unweighted_label
        vote_fraction[cell_index] = unweighted_score / k
        weighted_proportion[cell_index] = weighted_score / weights.sum()
        weighted_margin[cell_index] = (weighted_score - second_score) / weights.sum()

    return unweighted, weighted, vote_fraction, weighted_proportion, weighted_margin


def main():
    args = parse_args()
    if (args.output.exists() or args.summary.exists()) and not args.overwrite:
        raise FileExistsError("Output or summary exists. Use --overwrite to replace them.")
    temp_output = args.output.with_suffix(".tmp.h5ad")
    temp_summary = args.summary.with_suffix(".tmp.csv")

    target_names = np.load(args.work_dir / "target_obs_names.npy", allow_pickle=False)
    target_source_indices = np.load(
        args.work_dir / "target_source_indices.npy", allow_pickle=False
    )
    reference_names = np.load(
        args.work_dir / "reference_obs_names.npy", allow_pickle=False
    )
    reference_source_indices = np.load(
        args.work_dir / "reference_source_indices.npy", allow_pickle=False
    )
    if reference_source_indices.shape != (len(reference_names),):
        raise ValueError("Reference source-index mapping has the wrong shape.")
    if np.any(reference_source_indices < 0) or np.any(
        reference_source_indices >= args.expected_source_reference_cells
    ):
        raise ValueError("Reference source-index mapping is out of range.")
    target = ad.read_h5ad(args.target)
    if target.n_obs != args.expected_target_cells:
        raise ValueError(
            f"Expected {args.expected_target_cells} target cells, found {target.n_obs}."
        )
    if target_source_indices.shape != (len(target_names),):
        raise ValueError("Target source-index mapping has the wrong shape.")
    if np.any(target_source_indices < 0) or np.any(
        target_source_indices >= target.n_obs
    ):
        raise ValueError("Target source-index mapping is out of range.")
    if len(np.unique(target_source_indices)) != len(target_source_indices):
        raise ValueError("Target source-index mapping is not unique.")
    if not target.obs_names[target_source_indices].equals(pd.Index(target_names)):
        raise ValueError("Target AnnData order differs from prepared target features.")

    indices, distances, n_chunks = load_neighbors(
        args.chunks_dir,
        len(target_names),
        args.k,
    )
    if np.any(indices < 0) or np.any(indices >= len(reference_names)):
        raise ValueError("Combined neighbor indices are out of range.")
    if not np.isfinite(distances).all():
        raise ValueError("Combined neighbor distances contain non-finite values.")

    label_files = sorted(
        args.labels_dir.glob(f"res_*/{args.reference_section}.npy"),
        key=resolution_value,
    )
    if len(label_files) != EXPECTED_RESOLUTIONS:
        raise ValueError(
            f"Expected {EXPECTED_RESOLUTIONS} resolution label files, "
            f"found {len(label_files)}."
        )

    summary_rows = []
    obs_columns = {}
    n_full_target = target.n_obs
    for label_file in label_files:
        resolution_name = label_file.parent.name
        source_reference_labels = np.load(label_file, allow_pickle=False)
        if source_reference_labels.shape != (args.expected_source_reference_cells,):
            raise ValueError(f"Wrong reference label shape: {label_file}")
        reference_labels = source_reference_labels[reference_source_indices]

        results = transfer_one_resolution(reference_labels, indices, distances)
        unweighted, weighted, vote_fraction, weighted_proportion, weighted_margin = results
        categories = np.unique(
            np.concatenate(
                [reference_labels.astype(str), np.asarray([UNASSIGNED_LABEL], dtype=str)]
            )
        )

        full_unweighted = np.full(
            n_full_target, UNASSIGNED_LABEL, dtype=object
        )
        full_weighted = np.full(
            n_full_target, UNASSIGNED_LABEL, dtype=object
        )
        full_vote_fraction = np.full(n_full_target, np.nan, dtype=np.float32)
        full_weighted_proportion = np.full(n_full_target, np.nan, dtype=np.float32)
        full_weighted_margin = np.full(n_full_target, np.nan, dtype=np.float32)
        full_unweighted[target_source_indices] = unweighted.astype(str)
        full_weighted[target_source_indices] = weighted.astype(str)
        full_vote_fraction[target_source_indices] = vote_fraction
        full_weighted_proportion[target_source_indices] = weighted_proportion
        full_weighted_margin[target_source_indices] = weighted_margin

        obs_columns[f"transferred_{resolution_name}"] = pd.Categorical(
            full_unweighted, categories=categories
        )
        obs_columns[f"transferred_weighted_{resolution_name}"] = pd.Categorical(
            full_weighted, categories=categories
        )
        obs_columns[f"vote_fraction_{resolution_name}"] = full_vote_fraction
        obs_columns[f"weighted_vote_proportion_{resolution_name}"] = full_weighted_proportion
        obs_columns[f"weighted_top2_margin_{resolution_name}"] = full_weighted_margin

        disagreements = int(np.count_nonzero(unweighted != weighted))
        summary_rows.append(
            {
                "resolution": resolution_value(label_file),
                "resolution_name": resolution_name,
                "reference_n_types": len(np.unique(reference_labels)),
                "target_n_types_unweighted": len(np.unique(full_unweighted)),
                "target_n_types_weighted": len(np.unique(full_weighted)),
                "transferable_target_cells": len(unweighted),
                "unassigned_zero_shared_expression_cells": n_full_target - len(unweighted),
                "weighted_unweighted_disagreements": disagreements,
                "weighted_unweighted_disagreement_fraction": disagreements / len(unweighted),
                "mean_unweighted_vote_fraction": float(vote_fraction.mean()),
                "cells_unweighted_vote_below_0p5": int(
                    np.count_nonzero(vote_fraction < 0.5)
                ),
                "mean_weighted_vote_proportion": float(weighted_proportion.mean()),
                "mean_weighted_top2_margin": float(weighted_margin.mean()),
            }
        )
        print(
            f"transferred {resolution_name}: "
            f"{len(np.unique(reference_labels))} reference types, "
            f"{len(np.unique(unweighted))} target types"
        )

    full_indices = np.full((n_full_target, args.k), -1, dtype=np.int32)
    full_distances = np.full((n_full_target, args.k), np.nan, dtype=np.float32)
    full_indices[target_source_indices] = indices
    full_distances[target_source_indices] = distances
    obs_columns["knn_mean_distance"] = np.full(n_full_target, np.nan, dtype=np.float32)
    obs_columns["knn_min_distance"] = np.full(n_full_target, np.nan, dtype=np.float32)
    obs_columns["knn_mean_distance"][target_source_indices] = distances.mean(axis=1)
    obs_columns["knn_min_distance"][target_source_indices] = distances.min(axis=1)
    target.obs = pd.concat(
        [target.obs, pd.DataFrame(obs_columns, index=target.obs_names)],
        axis=1,
    )
    target.obsm["knn_reference_indices"] = full_indices
    target.obsm["knn_reference_distances"] = full_distances
    target.uns["label_transfer"] = {
        "k": int(args.k),
        "metric": "correlation",
        "algorithm": "brute",
        "primary_vote": "unweighted majority; ties resolved by mean neighbor distance",
        "normalization": "normalize_total(target_sum=10000), log1p",
        "resolutions": EXPECTED_RESOLUTIONS,
        "source_reference_cells": args.expected_source_reference_cells,
        "usable_reference_cells": len(reference_names),
        "excluded_zero_variance_reference_cells": int(
            args.expected_source_reference_cells - len(reference_names)
        ),
        "target_cells": int(n_full_target),
        "transferable_target_cells": int(len(target_names)),
        "excluded_zero_variance_target_cells": int(n_full_target - len(target_names)),
        "unassigned_label": UNASSIGNED_LABEL,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    temp_output.unlink(missing_ok=True)
    temp_summary.unlink(missing_ok=True)
    target.write_h5ad(temp_output)
    pd.DataFrame(summary_rows).sort_values("resolution").to_csv(
        temp_summary,
        index=False,
    )

    saved = ad.read_h5ad(temp_output, backed="r")
    try:
        if not saved.obs_names.equals(target.obs_names):
            raise ValueError("Saved target cell order changed.")
        if saved.obsm["knn_reference_indices"].shape != full_indices.shape:
            raise ValueError("Saved neighbor index shape is incorrect.")
    finally:
        saved.file.close()

    temp_output.replace(args.output)
    temp_summary.replace(args.summary)

    print(f"neighbor_chunks: {n_chunks}")
    print(f"output: {args.output}")
    print(f"summary: {args.summary}")


if __name__ == "__main__":
    main()
