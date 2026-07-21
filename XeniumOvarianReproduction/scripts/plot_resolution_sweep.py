#!/usr/bin/env python3
"""Plot raw Map Information across the ovarian Leiden resolution sweep."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "XeniumOvarianReproduction" / "results" / "full_map_info"
INPUT_CSV = RESULTS / "ovarian_percolation_scores_reduced.csv"
FIGURES_DIR = RESULTS / "figures"
RANKED_CSV = RESULTS / "raw_map_information_ranked.csv"

Y_COLUMN = "raw_score_weighted_mean"
BEST_COLUMNS = [
    "resolution",
    "total_number_of_types",
    "type_distribution_entropy",
    Y_COLUMN,
]


def load_results():
    """Load and validate the reduced resolution-sweep results."""
    data = pd.read_csv(INPUT_CSV)
    missing = [column for column in BEST_COLUMNS if column not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    for column in BEST_COLUMNS:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    invalid_rows = data[BEST_COLUMNS].isna().any(axis=1)
    if invalid_rows.any():
        raise ValueError(
            f"Found {int(invalid_rows.sum())} row(s) with missing or non-numeric "
            "required values"
        )
    if data.empty:
        raise ValueError(f"No results found in {INPUT_CSV}")
    if (data["resolution"] <= 0).any():
        raise ValueError("Resolution values must be positive for a logarithmic x-axis")

    return data


def plot_score_curve(data, x_column, x_label, title, output_name, log_x=False):
    """Plot one x-sorted raw Map Information curve and annotate its maximum."""
    sorted_data = data.sort_values(x_column, kind="stable")
    best = data.loc[data[Y_COLUMN].idxmax()]
    if log_x and (sorted_data[x_column] <= 0).any():
        raise ValueError(f"{x_column} values must be positive for a logarithmic x-axis")

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(
        sorted_data[x_column],
        sorted_data[Y_COLUMN],
        color="black",
        linewidth=2.0,
        marker="o",
        markersize=5,
    )
    if log_x:
        ax.set_xscale("log")
    ax.margins(y=0.18)

    maximum_label = (
        f"Global maximum\n"
        f"resolution: {best['resolution']:.6g}\n"
        f"total number of types: {best['total_number_of_types']:.6g}\n"
        f"type distribution entropy: {best['type_distribution_entropy']:.6g}\n"
        f"raw score weighted mean: {best[Y_COLUMN]:.6f}"
    )

    ax.scatter(
        [best[x_column]],
        [best[Y_COLUMN]],
        s=90,
        facecolors="white",
        edgecolors="black",
        linewidths=1.5,
        zorder=4,
    )
    ax.annotate(
        maximum_label,
        xy=(best[x_column], best[Y_COLUMN]),
        xytext=(0.98, 0.97),
        textcoords="axes fraction",
        fontsize=10,
        ha="right",
        va="top",
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "edgecolor": "0.7",
            "alpha": 0.95,
        },
    )
    ax.set_title(title, fontsize=16, pad=14)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel("Raw Map Information", fontsize=12)
    ax.grid(True, color="0.9", linewidth=0.8)
    fig.tight_layout()

    output = FIGURES_DIR / output_name
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output


def main():
    data = load_results()
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    output_paths = [
        plot_score_curve(
            data,
            x_column="resolution",
            x_label="Leiden resolution",
            title="Raw Map Information vs. Leiden resolution",
            output_name="raw_map_information_vs_resolution.png",
            log_x=True,
        ),
        plot_score_curve(
            data,
            x_column="total_number_of_types",
            x_label="Total number of types",
            title="Raw Map Information vs. number of types",
            output_name="raw_map_information_vs_number_of_types.png",
        ),
        plot_score_curve(
            data,
            x_column="total_number_of_types",
            x_label="Total number of types",
            title="Raw Map Information vs. number of types (log scale)",
            output_name="raw_map_information_vs_number_of_types_log.png",
            log_x=True,
        ),
        plot_score_curve(
            data,
            x_column="type_distribution_entropy",
            x_label="Type distribution entropy",
            title="Raw Map Information vs. type distribution entropy",
            output_name="raw_map_information_vs_type_entropy.png",
        ),
    ]

    ranked = data.sort_values(Y_COLUMN, ascending=False, kind="stable")
    ranked.to_csv(RANKED_CSV, index=False)
    best = ranked.iloc[0]

    print("Best row:")
    for column in BEST_COLUMNS:
        print(f"  {column}: {best[column]}")
    print("Generated outputs:")
    for path in [*output_paths, RANKED_CSV]:
        print(path)


if __name__ == "__main__":
    main()
