"""Fresh Xenium benchmark panels from saved native and comparison results.

Called by the matching R script, which supplies exact saved hexagon geometry.
Styles follow manuscript_update_20261012/generate_report.py; the map blocks
are arranged as the updated manuscript and Supplementary Figure S6.
"""
from pathlib import Path
import argparse
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="ccc_matplotlib_"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import FuncNorm, LinearSegmentedColormap, ListedColormap, Normalize, to_rgb
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from skimage.color import rgb2lab, lab2rgb

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path, required=True)
parser.add_argument("--geometry", type=Path, required=True)
args = parser.parse_args()
results = args.repo / "res/benchmarks/xenium"
figures = args.repo / "Figures"
METHODS = ["blisa", "stlearn", "spatialdm"]
LABELS = dict(zip(METHODS, ["blisa", "stLearn", "SpatialDM"]))
COLORS = dict(zip(METHODS, ["#222222", "#0072B2", "#D55E00"]))
REGIONS = ["Invasive_Tumor", "DCIS_1", "DCIS_2"]
REGION_LABELS = ["Invasive tumour", "DCIS 1", "DCIS 2"]
REGION_COLORS = ["#377EB8", "#4DAF4A", "#006D2C"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "pdf.fonttype": 42, "ps.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})

bins = pd.read_csv(args.geometry / "bins.csv")
vertices = pd.read_csv(args.geometry / "hexagons.csv")
borders = pd.read_csv(args.geometry / "roi_paths.csv")
hexagons = vertices[["x", "y"]].to_numpy().reshape(len(bins), 7, 2)
hex_limits = (hexagons[..., 0].min(), hexagons[..., 0].max(),
              hexagons[..., 1].min(), hexagons[..., 1].max())
summary = pd.read_csv(results / "xenium_benchmark_lr_summary.csv")
local = pd.read_csv(results / "xenium_benchmark_lri_by_bin.csv.gz", float_precision="round_trip", low_memory=False)
# The R wrapper converts the saved native object only for drawing these maps.
# Ignore legacy BLISA rows when reusing the existing Python-method cache.
blisa_local = pd.read_csv(args.geometry / "blisa_local.csv", float_precision="round_trip")
local = pd.concat([blisa_local, local[local.method != "blisa"]], ignore_index=True)
ids = summary.sort_values("rank").lr_id.drop_duplicates().tolist()
# 1. Read the shared LR set and regional proportions already computed in analysis.
overlap = pd.read_csv(results / "xenium_benchmark_jaccard_by_lr.csv")
common = set(overlap.loc[overlap.common_three_eligible, "lr_id"])
common_ids = [lr for lr in ids if lr in common]
regional = pd.read_csv(results / "xenium_benchmark_regional_coverage.csv", converters={"status": str})
region_lookup = regional.set_index(["lr_id", "region", "method"])
# Align the native local-result rows to the saved bin geometry for map drawing.
groups = {}
for (method, lr), group in local.groupby(["method", "lr_id"], sort=False):
    group = group.set_index("sample_id").loc[bins.sample_id]
    p = group.p_value.to_numpy(float)
    available = group.local_available.to_numpy(bool)
    selected = group.significant.to_numpy(bool)
    groups[(method, lr)] = (p, available, selected)
hot_cmap = LinearSegmentedColormap.from_list("hotspot", ["#FFFFCC", "#FFD700", "#FF7F00", "#D7301F"])


def save(fig, stem):
    fig.savefig(figures / f"{stem}.pdf", bbox_inches="tight", pad_inches=.04, facecolor="white")
    plt.close(fig)
    print(f"Rendered {stem}.pdf", flush=True)


def hex_map(ax, method, lr, header=False):
    p, available, selected = groups[(method, lr)]
    colors = np.tile(matplotlib.colors.to_rgba("#D0D0D0"), (len(bins), 1))
    colors[~available] = matplotlib.colors.to_rgba("#F2F2F2")
    colors[selected] = hot_cmap(Normalize(1, 4)(-np.log10(np.maximum(p[selected], 1e-4))))
    ax.add_collection(PolyCollection(hexagons, facecolors=colors, edgecolors="none",
        linewidths=0, antialiased=False, closed=True, rasterized=False))
    for (region, piece), border in borders.groupby(["region", "piece"], sort=False):
        ax.plot(border.x, border.y, color="white", lw=1.7)
        ax.plot(border.x, border.y, color=REGION_COLORS[REGIONS.index(region)], lw=1.1)
    if not available.any():
        ax.text(.5, .5, "NA", transform=ax.transAxes, ha="center", va="center", fontsize=9)
    ax.set_xlim(hex_limits[:2]); ax.set_ylim(hex_limits[2:]); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if header:
        ax.add_patch(plt.Rectangle((0, 1.015), 1, .145, transform=ax.transAxes,
            facecolor="#f2f3f5", edgecolor="none", clip_on=False))
        ax.plot([0, 1], [1.015, 1.015], transform=ax.transAxes,
            color=COLORS[method], lw=1.3, clip_on=False)
        ax.text(.5, 1.088, LABELS[method], transform=ax.transAxes,
            ha="center", va="center", fontsize=9, weight="bold")


def map_legends(fig, y, left=.14, width=.35, colorbar_left=.62, colorbar_width=.18):
    fig.legend(handles=[Patch(facecolor=c, label=label)
        for c, label in zip(REGION_COLORS, REGION_LABELS)],
        loc="center", bbox_to_anchor=(left + width / 2, y + .011),
        ncol=3, frameon=False, fontsize=8, handlelength=1.4, columnspacing=1.5)
    cax = fig.add_axes([colorbar_left, y, colorbar_width, .014])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(1, 4), cmap=hot_cmap),
        cax=cax, orientation="horizontal", ticks=[1, 2, 3, 4])
    cb.solids.set_rasterized(False)
    cb.ax.set_xticklabels(["1", "2", "3", "4+"])
    cb.ax.tick_params(labelsize=8, length=2)
    cb.outline.set_visible(False)
    fig.text(colorbar_left + colorbar_width + .025, y + .006,
        "−log₁₀\n(p-value)", ha="left", va="center", fontsize=8)


# 2. Figure 4g: selected LR, the three methods, and a shared horizontal legend.
fig = plt.figure(figsize=(7.5, 2.6))
grid = fig.add_gridspec(1, 3, left=.08, right=.995, top=.86, bottom=.20, wspace=.045)
axes = []
for j, method in enumerate(METHODS):
    ax = fig.add_subplot(grid[0, j])
    hex_map(ax, method, "CXCL12__CXCR4", header=True)
    axes.append(ax)
fig.canvas.draw()
box = axes[0].get_position()
fig.patches.append(plt.Rectangle((.025, box.y0), .028, box.height,
    transform=fig.transFigure, facecolor="#D0D0D0", edgecolor="none"))
fig.text(.039, box.y0 + box.height / 2, "CXCL12–CXCR4", rotation=90,
    ha="center", va="center", fontsize=9, style="italic", weight="bold")
map_legends(fig, .085, left=.03, width=.51, colorbar_left=.62, colorbar_width=.20)
save(fig, "figure4g")

# 3. Figure S6: all 12 common LR pairs, ordered as two blocks of six rows.
fig = plt.figure(figsize=(13, 9.2))
for block in range(2):
    left = .095 + block * .50
    grid = fig.add_gridspec(6, 3, left=left, right=left + .405,
        top=.95, bottom=.10, hspace=.08, wspace=.045)
    for i, lr in enumerate(common_ids[block * 6:(block + 1) * 6]):
        for j, method in enumerate(METHODS):
            ax = fig.add_subplot(grid[i, j])
            hex_map(ax, method, lr, header=i == 0)
            if j == 0:
                ax.text(-.11, .5, lr.replace("__", "–"), transform=ax.transAxes,
                    ha="right", va="center", fontsize=7.2, style="italic",
                    bbox=dict(facecolor="#e5e5e5", edgecolor="none", pad=3))
map_legends(fig, .06, left=.15, width=.35, colorbar_left=.59, colorbar_width=.17)
save(fig, "figureS6")

# 4. Figure S8a: untransformed percentages, shared log10(1 + percentage) colors.
coverage_rgb = np.array([to_rgb(c) for c in ["#FFFFFF", "#FFFFCC", "#FFD700", "#FF7F00", "#D7301F"]])
coverage_lab = rgb2lab(coverage_rgb.reshape(1, 5, 3)).reshape(5, 3)
lab_steps = np.column_stack([np.interp(np.linspace(0, 1, 1025), np.linspace(0, 1, 5),
    coverage_lab[:, i]) for i in range(3)])
coverage_cmap = ListedColormap(lab2rgb(lab_steps.reshape(1, -1, 3)).reshape(-1, 3))
coverage_cmap.set_bad("#D9D9D9")
coverage_norm = FuncNorm((lambda x: np.log10(1 + x), lambda y: 10 ** y - 1), vmin=0, vmax=100)
width, height = 6.44, 3.25
left, gap, bottom, cell_w, cell_h = 1.45, .16, .64, .40, .17
heat_h = len(common_ids) * cell_h
fig = plt.figure(figsize=(width, height))
for j, region in enumerate(REGIONS):
    total = int(region_lookup.loc[(common_ids[0], region, METHODS[0]), "n_total"])
    mat = np.array([[100 * region_lookup.loc[(lr, region, method), "proportion"]
        for method in METHODS] for lr in common_ids])
    ax = fig.add_axes([(left + j * (3 * cell_w + gap)) / width,
        bottom / height, 3 * cell_w / width, heat_h / height])
    mesh = ax.pcolormesh(np.arange(4) - .5, np.arange(len(common_ids) + 1) - .5,
        np.ma.masked_invalid(mat), norm=coverage_norm, cmap=coverage_cmap, shading="flat")
    ax.set_xlim(-.5, 2.5); ax.set_ylim(len(common_ids) - .5, -.5)
    ax.set_title(f"{REGION_LABELS[j]}\nn = {total:,}", fontsize=8,
        weight="bold", pad=8, color=REGION_COLORS[j])
    ax.set_xticks(range(3), [LABELS[m] for m in METHODS], rotation=50,
        ha="right", rotation_mode="anchor", fontsize=7)
    for tick, method in zip(ax.get_xticklabels(), METHODS):
        tick.set_color(COLORS[method]); tick.set_weight("bold")
    ax.set_yticks(range(len(common_ids)),
        [lr.replace("__", "–") for lr in common_ids] if j == 0 else [""] * len(common_ids),
        fontsize=7, weight="bold")
    for r in range(len(common_ids)):
        for c in range(3):
            value = mat[r, c]
            label = "/" if not np.isfinite(value) else ("<0.1%" if 0 < value < .1 else f"{value:.1f}%")
            luminance = np.dot(coverage_cmap(coverage_norm(value))[:3], [.2126, .7152, .0722]) if np.isfinite(value) else 1
            ax.text(c, r, label, ha="center", va="center", fontsize=6.4,
                color="white" if luminance < .48 else "#111111")
    ax.set_xticks(np.arange(-.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(common_ids), 1), minor=True)
    ax.grid(which="minor", color="#444444", linewidth=.45)
    ax.tick_params(which="both", length=0, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(True); spine.set_edgecolor("#111111"); spine.set_linewidth(.65)
cax = fig.add_axes([(width - .75) / width, (bottom + .1) / height,
    .10 / width, (heat_h - .2) / height])
ticks = [0, 1, 5, 10, 25, 50, 100]
cb = fig.colorbar(mesh, cax=cax, ticks=ticks)
cb.solids.set_rasterized(False); cb.solids.set_edgecolor("face")
cb.ax.set_yticklabels([f"{v}%" for v in ticks], fontsize=6.8)
cax.set_title("Region\ncoverage", fontsize=7, pad=8)
save(fig, "figureS8a")
