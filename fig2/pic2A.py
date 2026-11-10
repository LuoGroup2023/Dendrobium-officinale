#!/usr/bin/env python3
"""Plot the final collapsed-versus-haplotype DE comparison with matplotlib.

Figure contract
---------------
Core conclusion: tissue-effect directions are concordant between reference
representations, while the collapsed reference yields more significant calls.
Evidence: three matched-locus log2FC scatter plots, one per tissue contrast.
Archetype: quantitative grid.
Export: double-column SVG/PDF plus 600-dpi TIFF and PNG preview.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


CONTRASTS = (
    ("Stem_vs_Root", "Stem vs root"),
    ("Leaf_vs_Root", "Leaf vs root"),
    ("Stem_vs_Leaf", "Stem vs leaf"),
)

CATEGORY_ORDER = (
    "Neither significant",
    "Collapsed only",
    "Haplotype only",
    "Shared significant",
)

COLORS = {
    "Neither significant": "#B8BDC5",
    "Collapsed only": "#D08A45",
    "Haplotype only": "#2A8C82",
    "Shared significant": "#315F86",
}


def as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def category(row: dict[str, str]) -> str:
    collapsed = as_bool(row["collapsed_significant"])
    tetrad = as_bool(row["tetrad_significant"])
    if collapsed and tetrad:
        return "Shared significant"
    if collapsed:
        return "Collapsed only"
    if tetrad:
        return "Haplotype only"
    return "Neither significant"


def load_summary(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        return {row["contrast"]: row for row in csv.DictReader(handle, delimiter="\t")}


def load_detail(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if as_bool(row["tested_both"])
        ]


def write_source_data(path: Path, data: dict[str, list[dict[str, str]]]) -> None:
    columns = [
        "contrast",
        "TetradID",
        "collapsed_gene",
        "collapsed_log2FoldChange",
        "tetrad_log2FoldChange",
        "collapsed_padj",
        "tetrad_padj",
        "collapsed_significant",
        "tetrad_significant",
        "significance_category",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for contrast, rows in data.items():
            for row in rows:
                writer.writerow(
                    {
                        "contrast": contrast,
                        "TetradID": row["TetradID"],
                        "collapsed_gene": row["collapsed_gene"],
                        "collapsed_log2FoldChange": row["collapsed_log2FoldChange"],
                        "tetrad_log2FoldChange": row["tetrad_log2FoldChange"],
                        "collapsed_padj": row["collapsed_padj"],
                        "tetrad_padj": row["tetrad_padj"],
                        "collapsed_significant": row["collapsed_significant"],
                        "tetrad_significant": row["tetrad_significant"],
                        "significance_category": category(row),
                    }
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "data" / "collapsed_vs_haplotype_de",

    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent / "output" / "collapsed_vs_haplotype_de",
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    summary = load_summary(args.input_dir / "collapsed_vs_tetrad_DE_summary.tsv")
    data = {
        contrast: load_detail(
            args.input_dir / f"{contrast}__strict_one_to_one_comparison.tsv"
        )
        for contrast, _ in CONTRASTS
    }
    write_source_data(args.outdir / "Figure2C_source_data.tsv", data)

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "legend.frameon": False,
        }
    )

    width_in = 183 / 25.4
    height_in = 70 / 25.4
    fig, axes = plt.subplots(1, 3, figsize=(width_in, height_in), sharex=True, sharey=True)

    limit = 28.0
    ticks = [-20, -5, -1, 0, 1, 5, 20]
    scatter_sizes = {
        "Neither significant": 2.4,
        "Collapsed only": 4.2,
        "Haplotype only": 5.0,
        "Shared significant": 4.2,
    }
    alphas = {
        "Neither significant": 0.25,
        "Collapsed only": 0.58,
        "Haplotype only": 0.78,
        "Shared significant": 0.62,
    }

    for panel_index, (contrast, label) in enumerate(CONTRASTS):
        ax = axes[panel_index]
        rows = data[contrast]
        grouped = {name: [] for name in CATEGORY_ORDER}
        for row in rows:
            grouped[category(row)].append(row)

        ax.axhline(0, color="#D7D9DC", linewidth=0.55, zorder=0)
        ax.axvline(0, color="#D7D9DC", linewidth=0.55, zorder=0)
        ax.plot([-limit, limit], [-limit, limit], color="#555B61", linewidth=0.75, zorder=1)

        for name in CATEGORY_ORDER:
            subset = grouped[name]
            if not subset:
                continue
            x = np.fromiter(
                (float(row["collapsed_log2FoldChange"]) for row in subset), dtype=float
            )
            y = np.fromiter(
                (float(row["tetrad_log2FoldChange"]) for row in subset), dtype=float
            )
            ax.scatter(
                x,
                y,
                s=scatter_sizes[name],
                c=COLORS[name],
                alpha=alphas[name],
                linewidths=0,
                rasterized=True,
                zorder=2 if name == "Neither significant" else 3,
            )

        ax.set_xscale("symlog", linthresh=1, linscale=0.8, base=10)
        ax.set_yscale("symlog", linthresh=1, linscale=0.8, base=10)
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([str(value) for value in ticks])
        ax.set_yticklabels([str(value) for value in ticks])
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(label, pad=4, fontweight="bold")
        ax.set_xlabel("Collapsed-reference log$_2$FC")
        if panel_index == 0:
            ax.set_ylabel("Four-haplotype-summed log$_2$FC")

        stats = summary[contrast]
        annotation = (
            f"$n$ = {int(stats['tested_in_both']):,}\n"
            f"$r_s$ = {float(stats['log2fc_spearman']):.3f}\n"
            f"Shared DE = {int(stats['shared_significant']):,}\n"
            f"Opposite shared DE = {int(stats['shared_opposite_direction']):,}"
        )
        ax.text(
            0.035,
            0.965,
            annotation,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=6.1,
            linespacing=1.18,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 1.2},
        )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=COLORS[name],
            markeredgecolor="none",
            markersize=4,
            label=name,
        )
        for name in CATEGORY_ORDER
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=4,
        bbox_to_anchor=(0.5, 0.005),
        handletextpad=0.35,
        columnspacing=1.2,
        fontsize=6.3,
    )
    fig.subplots_adjust(left=0.075, right=0.995, top=0.93, bottom=0.22, wspace=0.25)

    prefix = args.outdir / "Figure2C_collapsed_vs_haplotype_DE"
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    legend_text = (

    )
    (args.outdir / "Figure2C_legend.txt").write_text(legend_text, encoding="ascii")

    with (args.outdir / "Figure2C_QA.txt").open("w", encoding="ascii") as handle:
        handle.write("backend\tPython matplotlib\n")
        handle.write("archetype\tquantitative grid\n")
        handle.write("width_mm\t183\n")
        handle.write("height_mm\t70\n")
        handle.write("axis_scale\tsymmetric log; linear threshold 1\n")
        handle.write("source_data\tFigure2C_source_data.tsv\n")
        for contrast, _ in CONTRASTS:
            handle.write(f"{contrast}_plotted_points\t{len(data[contrast])}\n")

    print("PASS: Figure 2C exports and source data created")
    print("DONE:", args.outdir)


if __name__ == "__main__":
    main()
