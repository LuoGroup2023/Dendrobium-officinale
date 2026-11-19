#!/usr/bin/env python3
"""Redraw Figure 2A-C with linear axes and accessible display styling.

Figure contract: compare matched tissue effects and observed versus technical
copy contributions without changing loci, numerical values, or statistics.
Archetype: quantitative grid, exported as separate replacement panels.
Python/matplotlib only; editable SVG/PDF, PNG preview and 600-dpi TIFF.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/doff_figure2_display_matplotlib")

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import numpy as np
import pandas as pd


REV = Path(__file__).resolve().parents[1]
ROOT = REV.parent
SAMPLES = tuple(f"{t}_R{r}" for t in ("Root", "Stem", "Leaf") for r in (1, 2, 3))
LABELS = tuple(f"{t}{r}" for t in ("Root", "Stem", "Leaf") for r in (1, 2, 3))
TISSUE_COLORS = {"Root": "#0072B2", "Stem": "#E69F00", "Leaf": "#CC79A7"}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(fig: plt.Figure, stem: Path) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        for label in ax.texts:
            box = label.get_window_extent(renderer)
            if box.x0 < 0 or box.y0 < 0 or box.x1 > fig.bbox.width or box.y1 > fig.bbox.height:
                raise ValueError(f"Annotation outside export canvas: {stem.name}: {label.get_text()}")
    for extension, kwargs in (
        ("svg", {}),
        ("pdf", {}),
        ("png", {"dpi": 300}),
        ("tiff", {"dpi": 600}),
    ):
        fig.savefig(stem.with_suffix("." + extension), facecolor="white", **kwargs)
    plt.close(fig)


def load_original_a():
    path = ROOT / "scripts/pic2a_road.py"
    spec = importlib.util.spec_from_file_location("original_figure2a", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def plot_a(out: Path, manifest: dict) -> None:
    original, original_path = load_original_a()
    source = ROOT / "work/collapsed_vs_haplotype_de"
    summary_path = source / "collapsed_vs_tetrad_DE_summary.tsv"
    summary = original.load_summary(summary_path)
    data = {
        contrast: original.load_detail(
            source / f"{contrast}__strict_one_to_one_comparison.tsv"
        )
        for contrast, _ in original.CONTRASTS
    }
    source_path = out / "Figure2A_source_data.tsv"
    original.write_source_data(source_path, data)
    previous = ROOT / "Figure2_final_sources/A_matched_log2FC/Figure2A_source_data.tsv"
    pd.testing.assert_frame_equal(pd.read_csv(source_path, sep="\t"),
                                  pd.read_csv(previous, sep="\t"))

    sizes = {"Neither significant": 2.4, "Collapsed only": 4.2,
             "Haplotype only": 5.0, "Shared significant": 4.2}
    alpha = {"Neither significant": .25, "Collapsed only": .58,
             "Haplotype only": .78, "Shared significant": .62}
    records = []
    with mpl.rc_context({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
                         "xtick.labelsize": 6, "ytick.labelsize": 6}):
        fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 70 / 25.4),
                                 sharex=True, sharey=True)
        fig.subplots_adjust(left=.095, right=.99, top=.91, bottom=.24, wspace=.28)
        for ax, (contrast, label) in zip(axes, original.CONTRASTS):
            rows = data[contrast]
            if len(rows) != 8592:
                raise ValueError(f"Unexpected Figure 2A population: {contrast}")
            ax.axhline(0, color="#D7D9DC", linewidth=.55, zorder=0)
            ax.axvline(0, color="#D7D9DC", linewidth=.55, zorder=0)
            ax.plot([-28, 28], [-28, 28], color="#555B61", linewidth=.75, zorder=1)
            for name in original.CATEGORY_ORDER:
                selected = [r for r in rows if original.category(r) == name]
                if not selected:
                    continue
                ax.scatter([float(r["collapsed_log2FoldChange"]) for r in selected],
                           [float(r["tetrad_log2FoldChange"]) for r in selected],
                           s=sizes[name], c=original.COLORS[name], alpha=alpha[name],
                           linewidths=0, rasterized=True,
                           zorder=2 if name == "Neither significant" else 3)
            ax.set_xscale("linear")
            ax.set_yscale("linear")
            ax.set_xlim(-28, 28)
            ax.set_ylim(-28, 28)
            ticks = [-20, -10, 0, 10, 20]
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.set_aspect("equal", adjustable="box")
            ax.set_title(label, pad=4, fontweight="bold")
            ax.set_xlabel("Collapsed-reference log$_2$FC")
            if ax is axes[0]:
                ax.set_ylabel("Four-haplotype-summed log$_2$FC")
            stats = summary[contrast]
            annotation = (f"$n$ = {int(stats['tested_in_both']):,}\n"
                          f"$r_s$ = {float(stats['log2fc_spearman']):.3f}\n"
                          f"Shared DE = {int(stats['shared_significant']):,}\n"
                          f"Opposite shared DE = {int(stats['shared_opposite_direction']):,}")
            ax.text(.035, .965, annotation, transform=ax.transAxes, ha="left", va="top",
                    fontsize=6.1, linespacing=1.18,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": .84, "pad": 1.2})
            values = np.array([[float(r[f]) for f in
                               ("collapsed_log2FoldChange", "tetrad_log2FoldChange")]
                              for r in rows])
            outside = int(np.any(np.abs(values) > 28, axis=1).sum())
            if outside:
                raise ValueError(f"Figure 2A would clip {outside} points in {contrast}")
            records.append({"contrast": contrast, "n": len(rows),
                            "x_min": float(values[:, 0].min()),
                            "x_max": float(values[:, 0].max()),
                            "y_min": float(values[:, 1].min()),
                            "y_max": float(values[:, 1].max()), "clipped_points": outside,
                            "x_scale": ax.get_xscale(), "y_scale": ax.get_yscale(),
                            "ticks": ticks})
        handles = [Line2D([0], [0], marker="o", linestyle="none",
                          markerfacecolor=original.COLORS[name], markeredgecolor="none",
                          markersize=4, label=name) for name in original.CATEGORY_ORDER]
        fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(.5, .005),
                   handletextpad=.35, columnspacing=1.2, fontsize=6.3, frameon=False)
        save(fig, out / "Figure2A_linear_axes")
    manifest["A"] = {"original_helper": str(original_path),
                     "original_helper_sha256": digest(original_path),
                     "input_sha256": {str(p): digest(p) for p in
                                      [summary_path] + [source / f"{c}__strict_one_to_one_comparison.tsv"
                                                        for c, _ in original.CONTRASTS]},
                     "source_table_matches_previous": True, "panels": records,
                     "changes": "Only symmetric-log to linear axes and equal integer major ticks"}


def plot_bc(out: Path, manifest: dict) -> None:
    path = REV / "results/20_uniform_count_statistics/aggregate/all_nine_technical_tests_with_global_FDR.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame = frame.loc[frame["testable"]].copy()
    if len(frame) != 2833 or frame.duplicated(["sample", "tetrad_id"]).any():
        raise ValueError("Unexpected Figure 2B/C population")
    if not np.allclose(frame["plugin_paired_excess"],
                       frame["observed_maximum_fraction"] - frame["plugin_null_median"]):
        raise ValueError("Paired differences do not reconstruct")
    columns = ["sample", "tetrad_id", "observed_maximum_fraction",
               "plugin_null_median", "plugin_paired_excess"]
    previous = REV / "figures/expression_revision"
    pd.testing.assert_frame_equal(frame[columns].reset_index(drop=True),
                                  pd.read_csv(previous / "Figure2BC_source_data.tsv", sep="\t"))
    frame[columns].to_csv(out / "Figure2BC_source_data.tsv", sep="\t", index=False)
    fig_b, ax_b = plt.subplots(figsize=(7.2, 3.0))
    fig_c, ax_c = plt.subplots(figsize=(7.2, 3.0))
    for fig in (fig_b, fig_c):
        fig.subplots_adjust(left=.12, right=.99, bottom=.19, top=.89)
    fig_b.subplots_adjust(top=.80)
    fig_c.subplots_adjust(top=.85)
    summary = []
    for x, (sample, label) in enumerate(zip(SAMPLES, LABELS)):
        color = TISSUE_COLORS[sample.split("_")[0]]
        chosen = frame.loc[frame["sample"] == sample]
        real = chosen["observed_maximum_fraction"].to_numpy()
        null = chosen["plugin_null_median"].to_numpy()
        excess = chosen["plugin_paired_excess"].to_numpy()
        for values, position, fill in ((null, x - .18, "#B9BDC1"), (real, x + .18, color)):
            violin = ax_b.violinplot(values, positions=[position], widths=.33,
                                     showextrema=False, points=160)
            for body in violin["bodies"]:
                body.set_facecolor(fill)
                body.set_edgecolor("none")
                body.set_alpha(.85)
            low, med, high = np.quantile(values, [.25, .5, .75])
            ax_b.plot([position, position], [low, high], color="#222222", lw=1.2, zorder=4)
            ax_b.scatter([position], [med], s=8, color="#111111", zorder=5)
        ax_b.text(x, 1.12, f"n={len(chosen)}", ha="center", va="bottom", fontsize=7.4,
                  transform=ax_b.get_xaxis_transform())
        ax_b.text(x, 1.035, f"{np.median(real):.3f}", ha="center", va="bottom",
                  fontsize=7.4, color="#222222", transform=ax_b.get_xaxis_transform())
        rng = np.random.default_rng(20261103 + x)
        ax_c.scatter(x + rng.uniform(-.22, .22, len(excess)), excess, s=5,
                     color=color, alpha=.4, linewidths=0, zorder=2)
        low, med, high = np.quantile(excess, [.25, .5, .75])
        ax_c.plot([x, x], [low, high], color="#111111", lw=2, zorder=4)
        ax_c.plot([x - .13, x + .13], [med, med], color="#111111", lw=1.5, zorder=5)
        ax_c.text(x, 1.025, f"{med:.3f}", ha="center", va="bottom", fontsize=7.8,
                  transform=ax_c.get_xaxis_transform())
        summary.append({"sample": sample, "n": len(chosen), "real_median": np.median(real),
                        "null_locus_median": np.median(null), "paired_excess_median": med,
                        "paired_excess_q25": low, "paired_excess_q75": high})
    for ax in (ax_b, ax_c):
        ax.set_xlim(-.6, 8.6)
        ax.set_xticks(range(9), LABELS)
        ax.tick_params(axis="both", direction="out")
        ax.set_facecolor("white")
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax_b.set_ylim(0, 1.01)
    ax_b.set_yticks([0, .25, .5, .75, 1.0])
    ax_b.set_ylabel("Maximum homolog fraction")
    values = frame["plugin_paired_excess"].to_numpy()
    low = min(-.05, float(values.min()) - .025)
    high = float(values.max()) + .035
    ax_c.set_ylim(low, high)
    ax_c.yaxis.set_major_locator(MultipleLocator(.10))
    ax_c.axhline(0, color="#555555", lw=.8, ls="--", zorder=1)
    ax_c.set_ylabel("Observed \u2212 null median", fontname="DejaVu Sans")
    ax_b.text(-.10, 1.08, "B", transform=ax_b.transAxes, fontsize=14, fontweight="bold")
    ax_c.text(-.10, 1.08, "C", transform=ax_c.transAxes, fontsize=14, fontweight="bold")
    new_summary = pd.DataFrame(summary)
    pd.testing.assert_frame_equal(new_summary,
                                  pd.read_csv(previous / "Figure2BC_summary.tsv", sep="\t"))
    new_summary.to_csv(out / "Figure2BC_summary.tsv", sep="\t", index=False)
    manifest["BC"] = {"source": str(path), "source_sha256": digest(path),
                      "source_table_matches_previous": True, "summary_matches_previous": True,
                      "n_locus_library_records": len(frame), "Primary_fixed": 534,
                      "minimum_assigned_fragments": 10,
                      "tissue_colors": TISSUE_COLORS, "axis_decimal_places": 2,
                      "B_limits": list(ax_b.get_ylim()), "C_limits": list(ax_c.get_ylim()),
                      "jitter_seed_unchanged": "20261103 + sample_index",
                      "changes": "Tick-label formatting, tissue colors, and dark annotation text only"}
    save(fig_b, out / "Figure2B_depth_matched_accessible")
    save(fig_c, out / "Figure2C_depth_matched_accessible")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path,
                        default=REV / "figures/Figure2_display_revision_20261112")
    args = parser.parse_args()
    if args.outdir.exists():
        raise FileExistsError(f"Use a new output directory: {args.outdir}")
    args.outdir.mkdir(parents=True)
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9,
                         "ytick.labelsize": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": .8, "axes.grid": False,
                         "svg.fonttype": "none", "pdf.fonttype": 42})
    manifest = {"purpose": "Display-only responses to suggestions 5-7",
                "backend": "Python/matplotlib", "new_statistical_tests": False,
                "script": str(Path(__file__).resolve()), "script_sha256": digest(Path(__file__)),
                "matplotlib": mpl.__version__, "numpy": np.__version__, "pandas": pd.__version__}
    plot_a(args.outdir, manifest)
    plot_bc(args.outdir, manifest)
    (args.outdir / "display_revision_audit.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("PASS: Figure 2A and B/C source tables match the previous panels")
    print("PASS: B/C counts, medians and IQRs unchanged")
    print("DONE:", args.outdir)


if __name__ == "__main__":
    main()
