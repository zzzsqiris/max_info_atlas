#!/usr/bin/env python3
"""Plot full ovarian Leiden cluster maps as bounding-box-clipped Voronoi cells."""

import argparse
import csv
import json
import resource
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from scipy.ndimage import binary_dilation
from scipy.spatial import Voronoi

from plot_full_map_info_results import (
    build_label_color_mapping,
    cluster_colors,
    format_resolution_dirname,
)


ROOT = Path(__file__).resolve().parents[2]
REPRODUCTION = ROOT / "XeniumOvarianReproduction"
DEFAULT_RESULTS_DIR = REPRODUCTION / "results" / "full_map_info"
XY_FILE = REPRODUCTION / "data" / "processed" / "xy_coordinates" / "ovarian_XY.npy"
DEFAULT_RESOLUTIONS = (0.095, 1.099)
MASK_PIXEL_SIZE = 10.0
MASK_RADIUS = 30.0


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot ovarian Leiden clusters as Voronoi cells."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Result directory containing clustering outputs and reduced scores",
    )
    parser.add_argument(
        "--resolutions",
        type=float,
        nargs="+",
        default=DEFAULT_RESOLUTIONS,
        help="One or more Leiden resolutions to plot",
    )
    return parser.parse_args()


def peak_rss_mib():
    """Return peak resident memory in MiB on macOS or Linux."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        peak *= 1024
    return peak / (1024**2)


def load_score_rows(scores_file, resolutions):
    """Load the requested raw scores and declared type counts from the sweep CSV."""
    if not scores_file.exists():
        raise FileNotFoundError(f"Reduced score CSV not found: {scores_file}")
    with scores_file.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    selected = {}
    for resolution in resolutions:
        matches = [
            row
            for row in rows
            if np.isclose(float(row["resolution"]), resolution, rtol=0, atol=5e-7)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one score row for resolution {resolution:.3f}, "
                f"found {len(matches)} in {scores_file}"
            )
        row = matches[0]
        selected[resolution] = {
            "raw_score": float(row["raw_score_weighted_mean"]),
            "number_of_types": int(float(row["total_number_of_types"])),
        }
    return selected


def load_inputs(clustering_dir, resolutions):
    """Load full XY and original labels without changing their row order."""
    xy = np.load(XY_FILE)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError(f"Expected XY shape (n_cells, 2), got {xy.shape}")
    if not np.isfinite(xy).all():
        raise ValueError(f"Non-finite coordinates found in {XY_FILE}")
    if len(np.unique(xy, axis=0)) != len(xy):
        raise ValueError("Duplicate XY coordinates cannot receive distinct Voronoi cells")

    labels_by_resolution = {}
    for resolution in resolutions:
        label_file = (
            clustering_dir
            / format_resolution_dirname(resolution)
            / "ovarian.npy"
        )
        if not label_file.exists():
            raise FileNotFoundError(
                f"Clustering labels not found for resolution {resolution}: {label_file}"
            )
        labels = np.load(label_file)
        if labels.ndim != 1:
            raise ValueError(f"Expected one-dimensional labels in {label_file}, got {labels.shape}")
        if len(labels) != len(xy):
            raise ValueError(
                f"Cell count mismatch: {label_file} has {len(labels):,} labels, "
                f"but XY has {len(xy):,} rows"
            )
        labels_by_resolution[resolution] = labels

    return xy, labels_by_resolution


def clip_polygon_to_rectangle(polygon, bounds):
    """Clip one ordered polygon to (xmin, xmax, ymin, ymax)."""
    xmin, xmax, ymin, ymax = bounds
    clipped = np.asarray(polygon, dtype=float)

    for axis, boundary, keep_above in (
        (0, xmin, True),
        (0, xmax, False),
        (1, ymin, True),
        (1, ymax, False),
    ):
        if len(clipped) == 0:
            break

        output = []
        previous = clipped[-1]
        previous_inside = previous[axis] >= boundary if keep_above else previous[axis] <= boundary

        for current in clipped:
            current_inside = current[axis] >= boundary if keep_above else current[axis] <= boundary
            if current_inside != previous_inside:
                fraction = (boundary - previous[axis]) / (current[axis] - previous[axis])
                intersection = previous + fraction * (current - previous)
                intersection[axis] = boundary
                output.append(intersection)
            if current_inside:
                output.append(current)
            previous = current
            previous_inside = current_inside

        clipped = np.asarray(output, dtype=float)

    return clipped


def finite_clipped_polygons(voronoi, bounds):
    """Return one finite, rectangle-clipped polygon per input point, in input order."""
    center = voronoi.points.mean(axis=0)
    radius = 2.0 * np.hypot(bounds[1] - bounds[0], bounds[3] - bounds[2])

    ridge_vertices = np.asarray(voronoi.ridge_vertices, dtype=int)
    infinite_ridge_mask = np.any(ridge_vertices < 0, axis=1)
    infinite_ridges = {}
    for points, vertices in zip(
        voronoi.ridge_points[infinite_ridge_mask],
        ridge_vertices[infinite_ridge_mask],
    ):
        p1, p2 = (int(points[0]), int(points[1]))
        v1, v2 = (int(vertices[0]), int(vertices[1]))
        if v1 >= 0:
            v1, v2 = v2, v1
        infinite_ridges.setdefault(p1, []).append((p2, v2))
        infinite_ridges.setdefault(p2, []).append((p1, v2))

    polygons = []
    for point_index, region_index in enumerate(voronoi.point_region):
        region = voronoi.regions[region_index]
        if not region:
            raise ValueError(f"Voronoi returned an empty region for cell {point_index}")

        finite_vertex_ids = [vertex for vertex in region if vertex >= 0]
        vertices = [voronoi.vertices[finite_vertex_ids]]

        if len(finite_vertex_ids) != len(region):
            far_vertices = []
            for neighbor_index, finite_vertex_id in infinite_ridges.get(point_index, []):
                tangent = voronoi.points[neighbor_index] - voronoi.points[point_index]
                tangent /= np.linalg.norm(tangent)
                normal = np.array([-tangent[1], tangent[0]])
                midpoint = (
                    voronoi.points[point_index] + voronoi.points[neighbor_index]
                ) / 2.0
                direction = normal if np.dot(midpoint - center, normal) >= 0 else -normal
                far_vertices.append(voronoi.vertices[finite_vertex_id] + direction * radius)
            if not far_vertices:
                raise ValueError(f"Could not close infinite Voronoi region for cell {point_index}")
            vertices.append(np.asarray(far_vertices))

        polygon = np.vstack(vertices)
        if len(finite_vertex_ids) != len(region):
            polygon_center = polygon.mean(axis=0)
            angles = np.arctan2(
                polygon[:, 1] - polygon_center[1],
                polygon[:, 0] - polygon_center[0],
            )
            polygon = polygon[np.argsort(angles)]

        if (
            (polygon[:, 0] < bounds[0]).any()
            or (polygon[:, 0] > bounds[1]).any()
            or (polygon[:, 1] < bounds[2]).any()
            or (polygon[:, 1] > bounds[3]).any()
        ):
            polygon = clip_polygon_to_rectangle(polygon, bounds)
        if len(polygon) < 3:
            raise ValueError(f"Clipping produced an invalid polygon for cell {point_index}")
        polygons.append(polygon)

    if len(polygons) != len(voronoi.points):
        raise AssertionError("Voronoi polygon order no longer matches input cell order")
    return polygons


def facecolors_for_labels(labels):
    """Map labels to colors without reordering cells."""
    unique_labels = np.sort(np.unique(labels))
    label_to_color = build_label_color_mapping(labels)
    return unique_labels, np.asarray([label_to_color[label] for label in labels])


def create_tissue_display_mask(xy, bounds):
    """Rasterize cells and dilate locally without filling enclosed tissue holes."""
    xmin, xmax, ymin, ymax = bounds
    width = max(1, int(np.ceil((xmax - xmin) / MASK_PIXEL_SIZE)))
    height = max(1, int(np.ceil((ymax - ymin) / MASK_PIXEL_SIZE)))

    x_index = np.floor((xy[:, 0] - xmin) / (xmax - xmin) * width).astype(int)
    y_index = np.floor((xy[:, 1] - ymin) / (ymax - ymin) * height).astype(int)
    x_index = np.clip(x_index, 0, width - 1)
    y_index = np.clip(y_index, 0, height - 1)

    occupied = np.zeros((height, width), dtype=bool)
    occupied[y_index, x_index] = True

    x_spacing = (xmax - xmin) / width
    y_spacing = (ymax - ymin) / height
    x_radius = int(np.ceil(MASK_RADIUS / x_spacing))
    y_radius = int(np.ceil(MASK_RADIUS / y_spacing))
    yy, xx = np.ogrid[-y_radius : y_radius + 1, -x_radius : x_radius + 1]
    footprint = (xx * x_spacing) ** 2 + (yy * y_spacing) ** 2 <= MASK_RADIUS**2
    return binary_dilation(occupied, structure=footprint)


def add_tissue_display_mask(ax, tissue_mask, bounds):
    """Hide polygons outside the tissue mask with a white rendering overlay."""
    white_overlay = np.ones((*tissue_mask.shape, 4), dtype=float)
    white_overlay[..., 3] = (~tissue_mask).astype(float)
    ax.imshow(
        white_overlay,
        extent=bounds,
        origin="lower",
        interpolation="nearest",
        aspect="auto",
        zorder=2,
    )


def plot_resolution(
    polygons, bounds, tissue_mask, labels, resolution, score_row, output_dir
):
    unique_labels, facecolors = facecolors_for_labels(labels)
    if len(unique_labels) != score_row["number_of_types"]:
        raise ValueError(
            f"Resolution {resolution:.3f} has {len(unique_labels)} labels, but "
            f"the score table declares {score_row['number_of_types']} types"
        )

    fig, ax = plt.subplots(figsize=(13, 8.5))
    collection = PolyCollection(
        polygons,
        facecolors=facecolors,
        edgecolors="none",
        linewidths=0,
        antialiaseds=False,
        rasterized=True,
    )
    ax.add_collection(collection)
    add_tissue_display_mask(ax, tissue_mask, bounds)
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X coordinate", fontsize=12)
    ax.set_ylabel("Y coordinate", fontsize=12)
    ax.set_title(
        f"Ovarian Leiden Voronoi clusters, resolution {resolution:.3f}\n"
        f"{len(unique_labels)} types | Raw Map Information: {score_row['raw_score']:.6f}",
        fontsize=16,
        pad=14,
    )
    _, counts = np.unique(labels, return_counts=True)
    handles = []
    for label, count, color in zip(unique_labels, counts, cluster_colors(len(unique_labels))):
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=color,
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

    output = output_dir / f"voronoi_clusters_{format_resolution_dirname(resolution)}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output


def record_stage(benchmark, name, started):
    elapsed = time.perf_counter() - started
    benchmark["stages"][name] = {
        "seconds": elapsed,
        "peak_rss_mib": peak_rss_mib(),
    }
    print(
        f"{name}: {elapsed:.2f} s; process peak RSS {benchmark['stages'][name]['peak_rss_mib']:.1f} MiB",
        flush=True,
    )


def main():
    args = parse_args()
    results_dir = args.results_dir
    resolutions = tuple(args.resolutions)
    clustering_dir = results_dir / "clustering" / "LeidenPCA50Correlation"
    scores_file = results_dir / "ovarian_percolation_scores_reduced.csv"
    output_dir = results_dir / "figures"
    benchmark_file = output_dir / "voronoi_benchmark.json"

    output_dir.mkdir(parents=True, exist_ok=True)
    total_started = time.perf_counter()
    benchmark = {
        "cell_count": None,
        "downsampled": False,
        "stages": {},
    }

    started = time.perf_counter()
    xy, labels_by_resolution = load_inputs(clustering_dir, resolutions)
    score_rows = load_score_rows(scores_file, resolutions)
    benchmark["cell_count"] = len(xy)
    record_stage(benchmark, "load_and_validate", started)

    bounds = (
        float(xy[:, 0].min()),
        float(xy[:, 0].max()),
        float(xy[:, 1].min()),
        float(xy[:, 1].max()),
    )

    started = time.perf_counter()
    voronoi = Voronoi(xy)
    record_stage(benchmark, "scipy_voronoi", started)

    started = time.perf_counter()
    polygons = finite_clipped_polygons(voronoi, bounds)
    benchmark["polygon_vertex_mib"] = sum(polygon.nbytes for polygon in polygons) / (1024**2)
    record_stage(benchmark, "reconstruct_and_clip", started)

    started = time.perf_counter()
    tissue_mask = create_tissue_display_mask(xy, bounds)
    benchmark["tissue_mask"] = {
        "method": "binary dilation of rasterized XY without hole filling",
        "pixel_size": MASK_PIXEL_SIZE,
        "radius": MASK_RADIUS,
        "shape": list(tissue_mask.shape),
        "tissue_fraction": float(tissue_mask.mean()),
    }
    record_stage(benchmark, "tissue_display_mask", started)

    outputs = []
    for resolution in resolutions:
        started = time.perf_counter()
        output = plot_resolution(
            polygons,
            bounds,
            tissue_mask,
            labels_by_resolution[resolution],
            resolution,
            score_rows[resolution],
            output_dir,
        )
        outputs.append(output)
        record_stage(benchmark, f"render_{resolution:.3f}", started)

    benchmark["total_seconds"] = time.perf_counter() - total_started
    benchmark["final_peak_rss_mib"] = peak_rss_mib()
    benchmark_file.write_text(json.dumps(benchmark, indent=2) + "\n")

    print("Generated outputs:")
    for output in outputs:
        print(output)
    print(benchmark_file)


if __name__ == "__main__":
    try:
        main()
    except MemoryError as error:
        raise SystemExit(
            "Full 407k-cell Voronoi exhausted memory; stopped without downsampling."
        ) from error
