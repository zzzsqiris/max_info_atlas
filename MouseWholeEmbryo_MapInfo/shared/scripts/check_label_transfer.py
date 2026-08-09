#!/usr/bin/env python3
"""Read-only completeness check for whole-embryo label-transfer chunks."""

import argparse
import json
from pathlib import Path

import numpy as np


PIPELINE_DIR = Path(__file__).resolve().parents[2]
VISIUM_RESULTS = PIPELINE_DIR / "datasets/visium_hd/results/label_transfer"
DEFAULT_WORK_DIR = VISIUM_RESULTS / "work"
DEFAULT_RESULTS_DIR = VISIUM_RESULTS / "chunk_results"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--k", type=int, default=15)
    return parser.parse_args()


def expected_ranges(path):
    ranges = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                start, stop = line.split("\t")
                ranges.append((int(start), int(stop)))
    return ranges


def main():
    args = parse_args()
    with open(args.work_dir / "manifest.json") as handle:
        manifest = json.load(handle)
    source_reference_cells = int(manifest["source_reference_cells"])
    reference = np.load(args.work_dir / "reference_features.npy", mmap_mode="r")
    reference_source_indices = np.load(
        args.work_dir / "reference_source_indices.npy", mmap_mode="r"
    )
    target = np.load(args.work_dir / "target_features.npy", mmap_mode="r")
    target_source_indices = np.load(
        args.work_dir / "target_source_indices.npy", mmap_mode="r"
    )
    if target_source_indices.shape != (len(target),):
        raise ValueError("Target source-index mapping has the wrong shape.")
    if reference_source_indices.shape != (len(reference),):
        raise ValueError("Reference source-index mapping has the wrong shape.")
    ranges = expected_ranges(args.work_dir / "jobs.txt")

    complete = 0
    missing = []
    invalid = []
    for start, stop in ranges:
        path = args.results_dir / f"neighbors_{start:06d}_{stop:06d}.npz"
        if not path.exists():
            missing.append(path.name)
            continue
        try:
            with np.load(path, allow_pickle=False) as chunk:
                shape = (stop - start, args.k)
                valid = (
                    int(chunk["start"]) == start
                    and int(chunk["stop"]) == stop
                    and chunk["indices"].shape == shape
                    and chunk["distances"].shape == shape
                    and np.all(chunk["indices"] >= 0)
                    and np.all(chunk["indices"] < len(reference))
                    and np.isfinite(chunk["distances"]).all()
                )
        except (OSError, ValueError, KeyError):
            valid = False
        if valid:
            complete += 1
        else:
            invalid.append(path.name)

    print(f"reference_cells: {len(reference)}")
    print(
        "excluded_zero_variance_reference_cells: "
        f"{source_reference_cells - len(reference)}"
    )
    full_target_cells = int(manifest["target_cells"])
    excluded_target_cells = int(manifest["excluded_zero_variance_target_cells"])
    if np.any(target_source_indices < 0) or np.any(target_source_indices >= full_target_cells):
        raise ValueError("Target source-index mapping is out of range.")
    if len(np.unique(target_source_indices)) != len(target_source_indices):
        raise ValueError("Target source-index mapping is not unique.")
    print(f"target_cells: {full_target_cells}")
    print(f"target_cells_for_knn: {len(target)}")
    print(f"excluded_zero_variance_target_cells: {excluded_target_cells}")
    print(f"expected_chunks: {len(ranges)}")
    print(f"complete_chunks: {complete}")
    print(f"missing_chunks: {len(missing)}")
    print(f"invalid_chunks: {len(invalid)}")
    print(f"progress: {100 * complete / len(ranges):.1f}%")
    for name in missing[:5]:
        print(f"  missing: {name}")
    for name in invalid[:5]:
        print(f"  invalid: {name}")

    if missing or invalid:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
