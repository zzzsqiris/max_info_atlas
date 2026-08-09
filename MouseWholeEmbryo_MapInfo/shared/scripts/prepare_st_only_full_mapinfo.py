#!/usr/bin/env python3
"""Prepare spatial inputs for the expression-only Map Information run."""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np

from max_info_atlases.percolation.graph_percolation import EdgeListManager


PIPELINE_DIR = Path(__file__).resolve().parents[2]
VISIUM_DIR = PIPELINE_DIR / "datasets/visium_hd"
DEFAULT_INPUT = VISIUM_DIR / "data/processed/mouse_embryo_visium_hd_raw.h5ad"
DEFAULT_INPUT_DIR = VISIUM_DIR / "data/processed/visium_hd_st_only_mapinfo"
DEFAULT_EDGE_DIR = VISIUM_DIR / "results/st_only_map_info/edge_lists"
SECTION = "MouseEmbryo_VisiumHD"
EXPECTED_CELLS = 449_700


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--edge-list-dir", type=Path, default=DEFAULT_EDGE_DIR)
    parser.add_argument("--section", default=SECTION)
    parser.add_argument("--expected-cells", type=int, default=EXPECTED_CELLS)
    parser.add_argument("--max-k", type=int, default=500)
    parser.add_argument(
        "--build-edge-list",
        action="store_true",
        help="Build the spatial edge list in addition to exporting Sections/XY.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def save_or_validate(path, values, overwrite):
    if path.exists() and not overwrite:
        saved = np.load(path, allow_pickle=False)
        if np.array_equal(saved, values):
            print(f"already valid: {path}")
            return
        raise ValueError(f"Existing file differs: {path}. Use --overwrite.")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, values, allow_pickle=False)
    print(f"saved: {path}")


def main():
    args = parse_args()
    if args.max_k < 1:
        raise ValueError("--max-k must be positive.")

    adata = ad.read_h5ad(args.input, backed="r")
    try:
        if adata.n_obs != args.expected_cells:
            raise ValueError(f"Expected {args.expected_cells} cells, found {adata.n_obs}.")
        if "spatial" not in adata.obsm:
            raise ValueError("AnnData is missing obsm['spatial'].")
        if "sample" not in adata.obs:
            raise ValueError("AnnData is missing obs['sample'].")
        sections = adata.obs["sample"].astype(str).to_numpy(dtype=str)
        xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    finally:
        adata.file.close()

    if set(sections) != {args.section}:
        raise ValueError(
            f"Expected one section named {args.section}, found {set(sections)}."
        )
    if xy.shape != (args.expected_cells, 2) or not np.isfinite(xy).all():
        raise ValueError(f"Invalid spatial coordinates: {xy.shape}.")
    if len(np.unique(xy, axis=0)) != len(xy):
        raise ValueError("Spatial coordinates contain duplicates.")

    save_or_validate(args.output_dir / "Sections.npy", sections, args.overwrite)
    save_or_validate(
        args.output_dir / "xy_coordinates" / f"{args.section}_XY.npy", xy, args.overwrite
    )

    if args.build_edge_list:
        edge_path = args.edge_list_dir / f"{args.section}_k{args.max_k}.npz"
        if edge_path.exists() and not args.overwrite:
            print(f"edge list already exists: {edge_path}")
        else:
            args.edge_list_dir.mkdir(parents=True, exist_ok=True)
            manager = EdgeListManager(base_dir=str(args.edge_list_dir))
            edges = manager.compute_edge_list(xy, args.max_k)
            manager.save_edge_list(xy, args.max_k, args.section, edges)
            print(f"saved edge list: {edge_path}")
            print(f"edge_list_shape: {edges.shape}")

    print(f"cells: {len(xy)}")
    print(f"section: {args.section}")


if __name__ == "__main__":
    main()
