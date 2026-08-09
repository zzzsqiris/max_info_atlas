#!/usr/bin/env python3
"""Compute exact correlation kNN for spatial target-cell ranges."""

import argparse
from pathlib import Path

import numpy as np
from sklearn.neighbors import NearestNeighbors


PIPELINE_DIR = Path(__file__).resolve().parents[2]
DEFAULT_WORK_DIR = PIPELINE_DIR / "datasets/visium_hd/results/label_transfer/work"
DEFAULT_OUTPUT_DIR = (
    PIPELINE_DIR / "datasets/visium_hd/results/label_transfer/chunk_results"
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-file", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--k", type=int, default=15)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def read_ranges(path):
    ranges = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(f"Expected start and stop in {path}: {line}")
            ranges.append((int(fields[0]), int(fields[1])))
    if not ranges:
        raise ValueError(f"No target ranges found in {path}.")
    return ranges


def valid_existing_output(path, start, stop, k, n_reference):
    """Return True only when an existing chunk is complete and usable."""
    try:
        with np.load(path, allow_pickle=False) as chunk:
            saved_start = int(chunk["start"])
            saved_stop = int(chunk["stop"])
            indices = chunk["indices"]
            distances = chunk["distances"]
    except (OSError, ValueError, KeyError):
        return False

    expected_shape = (stop - start, k)
    return (
        saved_start == start
        and saved_stop == stop
        and indices.shape == expected_shape
        and distances.shape == expected_shape
        and np.all(indices >= 0)
        and np.all(indices < n_reference)
        and np.isfinite(distances).all()
    )


def main():
    args = parse_args()
    if args.k < 1:
        raise ValueError("k must be positive.")

    reference = np.load(args.work_dir / "reference_features.npy", mmap_mode="r")
    target = np.load(args.work_dir / "target_features.npy", mmap_mode="r")
    if args.k > len(reference):
        raise ValueError("k exceeds the reference cell count.")

    neighbors = NearestNeighbors(
        n_neighbors=args.k,
        metric="correlation",
        algorithm="brute",
        n_jobs=1,
    )
    neighbors.fit(reference)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for start, stop in read_ranges(args.chunk_file):
        if start < 0 or stop <= start or stop > len(target):
            raise ValueError(f"Invalid target range: {start}:{stop}.")
        output = args.output_dir / f"neighbors_{start:06d}_{stop:06d}.npz"
        if output.exists() and not args.overwrite:
            if valid_existing_output(output, start, stop, args.k, len(reference)):
                print(f"already complete: {output}")
                continue
            raise ValueError(
                f"Existing output is incomplete or invalid: {output}. "
                "Use --overwrite to replace it."
            )

        distances, indices = neighbors.kneighbors(
            target[start:stop],
            return_distance=True,
        )
        distances = np.maximum(distances, 0.0)
        if distances.shape != (stop - start, args.k):
            raise ValueError("Neighbor distance shape is incorrect.")
        if not np.isfinite(distances).all():
            raise ValueError("Neighbor distances contain non-finite values.")
        if np.any(indices < 0) or np.any(indices >= len(reference)):
            raise ValueError("Neighbor indices are out of range.")

        np.savez_compressed(
            output,
            start=np.int64(start),
            stop=np.int64(stop),
            indices=indices.astype(np.int32, copy=False),
            distances=distances.astype(np.float32, copy=False),
        )
        print(f"saved: {output}")


if __name__ == "__main__":
    main()
