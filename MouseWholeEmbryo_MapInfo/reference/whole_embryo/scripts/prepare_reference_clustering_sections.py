#!/usr/bin/env python3
"""Create the one-section array required for Qiu E15.5 reference clustering."""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np


TASK_DIR = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    TASK_DIR / "data/processed/qiu_e15_5_whole_embryo_reference_raw.h5ad"
)
DEFAULT_OUTPUT = (
    TASK_DIR / "data/processed/qiu_e15_5_whole_embryo_reference_sections.npy"
)
EXPECTED_CELLS = 215_284
EXPECTED_DAY = "E15.5"
EXPECTED_EMBRYO = "embryo_59"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(
            f"Output already exists: {args.output}. Use --overwrite to replace it."
        )

    adata = ad.read_h5ad(args.input, backed="r")
    try:
        if adata.n_obs != EXPECTED_CELLS:
            raise ValueError(f"Expected {EXPECTED_CELLS} cells, found {adata.n_obs}.")
        for column in ["sample", "day", "embryo_id"]:
            if column not in adata.obs.columns:
                raise ValueError(f"AnnData obs does not contain '{column}'.")
        if set(adata.obs["day"].astype(str)) != {EXPECTED_DAY}:
            raise ValueError("Reference contains an unexpected developmental day.")
        if set(adata.obs["embryo_id"].astype(str)) != {EXPECTED_EMBRYO}:
            raise ValueError("Reference contains an unexpected embryo.")

        sections = adata.obs["sample"].astype(str).to_numpy(dtype=str)
        unique_sections = np.unique(sections)
        if len(unique_sections) != 1:
            raise ValueError(f"Expected one sample, found: {unique_sections.tolist()}.")
    finally:
        adata.file.close()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.output, sections, allow_pickle=False)

    saved = np.load(args.output, allow_pickle=False)
    if not np.array_equal(saved, sections):
        raise ValueError("Saved section array does not match AnnData cell order.")

    print(f"output: {args.output}")
    print(f"cells: {len(saved)}")
    print(f"section: {unique_sections[0]}")


if __name__ == "__main__":
    main()
