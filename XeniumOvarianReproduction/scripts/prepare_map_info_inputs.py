#!/usr/bin/env python3
"""Prepare the spatial input files expected by the professor's pipeline.

This script does not normalize expression, calculate PCA, or modify the source
AnnData file. It only exports:

    data/processed/Sections.npy
    data/processed/xy_coordinates/<section>_XY.npy

Optionally, ``--small-n-cells`` creates a spatially contiguous AnnData subset
for a local smoke test. That subset is only for checking that the pipeline can
run; it must not be used as the final full-data reproduction result.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


REPRODUCTION_DIR = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = REPRODUCTION_DIR / "data" / "raw" / "ovarian_xenium_genes.h5ad"
DEFAULT_OUTPUT = REPRODUCTION_DIR / "data" / "processed"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export section labels and XY coordinates for Map Info."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Source .h5ad file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output directory (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--section-name",
        default="ovarian",
        help="Section label assigned to every cell (default: ovarian)",
    )
    parser.add_argument(
        "--small-n-cells",
        type=int,
        default=None,
        help=(
            "Also create a spatially contiguous subset .h5ad with this many "
            "cells for a local smoke test"
        ),
    )
    return parser.parse_args()


def read_spatial(input_path: Path) -> np.ndarray:
    """Read obsm['spatial'] without loading the expression matrix."""
    with h5py.File(input_path, "r") as handle:
        if "obsm" not in handle or "spatial" not in handle["obsm"]:
            raise KeyError("The .h5ad file does not contain adata.obsm['spatial']")
        spatial = np.asarray(handle["obsm/spatial"])

    if spatial.ndim != 2 or spatial.shape[1] != 2:
        raise ValueError(f"Expected spatial shape (n_cells, 2), got {spatial.shape}")
    if not np.isfinite(spatial).all():
        raise ValueError("Spatial coordinates contain NaN or infinite values")
    return spatial


def save_pipeline_inputs(
    spatial: np.ndarray,
    output_dir: Path,
    section_name: str,
) -> None:
    xy_dir = output_dir / "xy_coordinates"
    xy_dir.mkdir(parents=True, exist_ok=True)

    sections = np.full(
        spatial.shape[0], section_name, dtype=f"<U{len(section_name)}"
    )
    np.save(output_dir / "Sections.npy", sections)
    np.save(xy_dir / f"{section_name}_XY.npy", spatial)


def create_small_subset(
    input_path: Path,
    output_dir: Path,
    section_name: str,
    spatial: np.ndarray,
    n_cells: int,
) -> None:
    if n_cells <= 0:
        raise ValueError("--small-n-cells must be greater than zero")
    if n_cells >= len(spatial):
        raise ValueError(
            f"--small-n-cells must be smaller than the full dataset ({len(spatial)})"
        )

    try:
        import anndata as ad
    except ImportError as exc:
        raise RuntimeError(
            "Creating a small subset requires anndata. Activate the professor's "
            "max_info_atlases environment first."
        ) from exc

    # Choose cells nearest the spatial median so the test set is one local region.
    center = np.median(spatial, axis=0)
    squared_distance = np.sum((spatial - center) ** 2, axis=1)
    selected = np.argpartition(squared_distance, n_cells - 1)[:n_cells]
    selected.sort()  # Preserve the source AnnData row order.

    print(f"Loading AnnData to create a {n_cells:,}-cell local test subset ...")
    adata = ad.read_h5ad(input_path)
    subset = adata[selected].copy()
    subset.obs["section"] = section_name

    small_dir = output_dir / f"small_{n_cells}"
    small_dir.mkdir(parents=True, exist_ok=True)
    subset_path = small_dir / f"ovarian_small_{n_cells}.h5ad"
    subset.write_h5ad(subset_path)

    subset_spatial = np.asarray(subset.obsm["spatial"])
    save_pipeline_inputs(subset_spatial, small_dir, section_name)

    print(f"Small test AnnData: {subset_path}")
    print(f"Small test Sections: {small_dir / 'Sections.npy'}")
    print(f"Small test XY: {small_dir / 'xy_coordinates' / f'{section_name}_XY.npy'}")


def main() -> None:
    args = parse_args()
    input_path = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not input_path.is_file():
        raise FileNotFoundError(f"Input .h5ad not found: {input_path}")
    if not args.section_name or "/" in args.section_name:
        raise ValueError("--section-name must be a non-empty filename-safe label")

    spatial = read_spatial(input_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    save_pipeline_inputs(spatial, output_dir, args.section_name)

    print(f"Source AnnData (unchanged): {input_path}")
    print(f"Cells: {len(spatial):,}")
    print(f"Sections: {output_dir / 'Sections.npy'}")
    print(
        "XY: "
        f"{output_dir / 'xy_coordinates' / f'{args.section_name}_XY.npy'}"
    )

    if args.small_n_cells is not None:
        create_small_subset(
            input_path=input_path,
            output_dir=output_dir,
            section_name=args.section_name,
            spatial=spatial,
            n_cells=args.small_n_cells,
        )


if __name__ == "__main__":
    main()
