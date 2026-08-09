#!/usr/bin/env python3
"""Create full-data figures that directly support Map Information."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = (
    ROOT / "XeniumOvarianReproduction" / "results" / "full_map_info"
)
XY_FILE = (
    ROOT
    / "XeniumOvarianReproduction"
    / "data"
    / "processed"
    / "xy_coordinates"
    / "ovarian_XY.npy"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot ovarian Map Information results."
    )
    parser.add_argument(
        "--resolution",
        required=True,
        type=float,
        help="Leiden resolution to use for the spatial cluster plot",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Result directory containing clustering and percolation outputs",
    )
    return parser.parse_args()


def format_resolution_dirname(resolution):
    return f"res_{resolution:.3f}".replace(".", "p")


def load_results(percolation_file):
    with np.load(percolation_file) as data:
        return {
            "xy": data["XY"],
            "labels": data["type_vec"],
            "ent_real": data["ent_real"],
            "ent_perm": data["ent_perm"],
            "pbond": data["pbond_vec"],
            "max_k": int(data["maxK"]),
        }


def load_scores(score_file):
    raw, normalized = score_file.read_text().strip().split("\t")
    return float(raw), float(normalized)


def load_spatial_results(resolution, clustering_dir):
    resolution_dirname = format_resolution_dirname(resolution)
    labels_file = clustering_dir / resolution_dirname / "ovarian.npy"

    if not XY_FILE.exists():
        raise FileNotFoundError(f"Spatial XY file not found: {XY_FILE}")
    if not labels_file.exists():
        raise FileNotFoundError(
            f"Clustering labels not found for resolution {resolution}: {labels_file}"
        )

    xy = np.load(XY_FILE)
    labels = np.load(labels_file)

    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError(f"Expected XY shape (N, 2), got {xy.shape} from {XY_FILE}")
    if labels.ndim != 1:
        raise ValueError(
            f"Expected one-dimensional clustering labels, got {labels.shape} "
            f"from {labels_file}"
        )
    if xy.shape[0] != labels.shape[0]:
        raise ValueError(
            f"XY/label length mismatch: XY has {xy.shape[0]} rows, "
            f"labels has {labels.shape[0]} values"
        )

    return {
        "xy": xy,
        "labels": labels,
        "resolution": resolution,
        "resolution_dirname": resolution_dirname,
    }


def plot_entropy_curves(
    results,
    raw_score,
    normalized_score,
    output_dir,
    title="Full ovarian dataset: percolation entropy curves",
    output_stem="full_percolation_entropy_curves",
):
    ent_real = results["ent_real"]
    ent_perm = results["ent_perm"]
    threshold = results["pbond"][: len(ent_real)]

    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    calculated_raw = float(
        integrate(np.abs(ent_perm - ent_real), x=threshold)
    )
    if not np.isclose(calculated_raw, raw_score, rtol=1e-9, atol=1e-9):
        raise ValueError(
            f"Curve integral {calculated_raw} does not match score {raw_score}"
        )

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(
        threshold,
        ent_real,
        color="black",
        linewidth=2.2,
        label="Real cluster labels",
    )
    ax.plot(
        threshold,
        ent_perm,
        color="0.45",
        linewidth=2.2,
        linestyle="--",
        label="Permuted cluster labels",
    )
    ax.fill_between(
        threshold,
        ent_real,
        ent_perm,
        color="0.85",
        alpha=0.8,
        label="Absolute difference integrated for raw score",
    )

    ax.set_title(title, fontsize=16, pad=14)
    ax.set_xlabel("Edge activation threshold", fontsize=12)
    ax.set_ylabel("Connected-component entropy (bits)", fontsize=12)
    ax.set_xlim(threshold.min(), threshold.max())
    ax.grid(True, color="0.9", linewidth=0.8)
    ax.legend(frameon=False, fontsize=10, loc="best")
    ax.text(
        0.02,
        0.04,
        (
            f"Cells: {len(results['labels']):,}\n"
            f"maxK: {results['max_k']}\n"
            f"Raw score: {raw_score:.6f}\n"
            f"Normalized score: {normalized_score:.6f}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
    )
    fig.tight_layout()
    outputs = [output_dir / f"{output_stem}.png", output_dir / f"{output_stem}.pdf"]
    for output in outputs:
        fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return outputs


def cluster_colors(n_clusters):
    tab20 = list(plt.get_cmap("tab20").colors)
    tab20b = list(plt.get_cmap("tab20b").colors)
    colors = tab20 + tab20b
    if n_clusters > len(colors):
        raise ValueError(
            f"Deterministic palette supports at most {len(colors)} clusters; "
            f"found {n_clusters}"
        )
    return colors[:n_clusters]


def build_label_color_mapping(labels):
    """Map sorted cluster labels to deterministic colors."""
    unique_labels = np.sort(np.unique(labels))
    colors = cluster_colors(len(unique_labels))
    return {label: color for label, color in zip(unique_labels, colors)}


def plot_spatial_clusters(results, output_dir):
    xy = results["xy"]
    labels = results["labels"]
    resolution = results["resolution"]
    resolution_dirname = results["resolution_dirname"]
    unique_labels, counts = np.unique(labels, return_counts=True)
    label_to_color = build_label_color_mapping(labels)
    point_colors = np.asarray([label_to_color[label] for label in labels])

    # Shuffle only the plotting order so one cluster is not systematically drawn last.
    order = np.random.default_rng(42).permutation(len(labels))

    fig, ax = plt.subplots(figsize=(13, 8.5))
    ax.scatter(
        xy[order, 0],
        xy[order, 1],
        c=point_colors[order],
        s=0.22,
        alpha=0.75,
        linewidths=0,
        rasterized=True,
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X coordinate", fontsize=12)
    ax.set_ylabel("Y coordinate", fontsize=12)
    ax.set_title(
        "Full ovarian dataset: Leiden spatial clusters\n"
        f"Resolution: {resolution:.3f} | Clusters: {len(unique_labels):,} | "
        f"Cells: {len(labels):,}",
        fontsize=16,
        pad=14,
    )
    handles = []
    for label, count in zip(unique_labels, counts):
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=label_to_color[label],
                markeredgecolor="none",
                markersize=6,
                label=f"Cluster {label} ({count:,})",
            )
        )
    ax.legend(
        handles=handles,
        title="Leiden cluster (cell count)",
        frameon=False,
        fontsize=8,
        title_fontsize=9,
        ncol=2,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.tight_layout()
    output = output_dir / f"full_spatial_leiden_clusters_{resolution_dirname}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output


def main():
    args = parse_args()
    results_dir = args.results_dir
    output_dir = results_dir / "figures"
    clustering_dir = results_dir / "clustering" / "LeidenPCA50Correlation"
    percolation_file = (
        results_dir
        / "percolation_results"
        / "LeidenPCA50Correlation"
        / "res_1p000"
        / "ovarian.npy.npz"
    )
    score_file = percolation_file.with_name("ovarian.score")

    output_dir.mkdir(parents=True, exist_ok=True)
    spatial_results = load_spatial_results(args.resolution, clustering_dir)
    spatial_path = plot_spatial_clusters(spatial_results, output_dir)

    print(f"XY shape: {spatial_results['xy'].shape}")
    print(f"Label shape: {spatial_results['labels'].shape}")
    print(f"Clusters: {len(np.unique(spatial_results['labels']))}")
    print(spatial_path)

    if percolation_file.exists() and score_file.exists():
        results = load_results(percolation_file)
        raw_score, normalized_score = load_scores(score_file)
        entropy_paths = plot_entropy_curves(
            results, raw_score, normalized_score, output_dir
        )
        for entropy_path in entropy_paths:
            print(entropy_path)
    else:
        print(
            "Skipping percolation entropy plot because its existing inputs are missing: "
            f"{percolation_file}, {score_file}"
        )


if __name__ == "__main__":
    main()
