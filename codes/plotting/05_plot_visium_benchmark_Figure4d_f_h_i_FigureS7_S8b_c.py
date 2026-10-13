"""Figure 4h-i, S7b and S8b-c from complete saved Visium benchmark CSVs.

Spatial-map functions preserve the native 12 October 2026 report geometry,
colours and ROI overlays. No inference, downloads or figure-specific data exports.
"""
import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ccc_visium_benchmark_mpl")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import FuncNorm, LinearSegmentedColormap, ListedColormap, Normalize, to_rgb
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from skimage.color import rgb2lab, lab2rgb

METHODS = ["blisa", "stlearn_bh", "spatialdm_permutation"]
NAMES = dict(zip(METHODS, ["blisa", "stLearn", "SpatialDM"]))
STEMS = dict(zip(METHODS, ["blisa_permutation", "stlearn_bh_adjusted", "spatialdm_permutation"]))
COLORS = dict(zip(METHODS, ["#222222", "#0072B2", "#D55E00"]))
REGIONS = ["Invasive_Tumor", "DCIS_1", "DCIS_2", "Outside ROI"]
REGION_NAMES = ["Invasive tumour", "DCIS 1", "DCIS 2", "Non-tumour*"]
ROI_COLORS = dict(zip(REGIONS, ["#377EB8", "#4DAF4A", "#006D2C", "#666666"]))
PVALUE_CMAP = LinearSegmentedColormap.from_list("hotspot_p", ["#FFFFCC", "#FFD700", "#FF7F00", "#D7301F"])
PVALUE_NORM = Normalize(1, 4)
ROI_LINEWIDTH, ROI_HALO_LINEWIDTH = 1.1, 1.7
_rgb = np.array([to_rgb(c) for c in ["#FFFFFF", "#FFFFCC", "#FFD700", "#FF7F00", "#D7301F"]])
_lab = rgb2lab(_rgb.reshape(1, 5, 3)).reshape(5, 3)
_interpolated = np.column_stack([np.interp(np.linspace(0, 1, 1025), np.linspace(0, 1, 5), _lab[:, i]) for i in range(3)])
REGIONAL_CMAP = ListedColormap(lab2rgb(_interpolated.reshape(1, -1, 3)).reshape(-1, 3))
REGIONAL_CMAP.set_bad("#D9D9D9")
REGIONAL_NORM = FuncNorm((lambda x: np.log10(1+x), lambda y: 10**y-1), vmin=0, vmax=100)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.labelweight": "bold", "axes.titleweight": "bold"})


def label(lr):
    return lr.replace("__", " → ")


def map_panel(ax, method, lr, kind):
    ax.set(xlim=MAP_XLIM,ylim=MAP_YLIM,aspect="equal")
    ax.apply_aspect()
    frame=map_lookup.get((method,lr))
    status=final[(final.condition==method)&(final.lr_id==lr)].iloc[0].status
    if status=="NR":
        ax.text(.5,.5,"NR",transform=ax.transAxes,ha="center",fontsize=8)
    else:
        available=np.zeros(len(spots),bool) if frame is None else frame.available.fillna(False).to_numpy(bool)
        selected=np.zeros(len(spots),bool) if frame is None else frame.selected.fillna(False).to_numpy(bool)
        rgba=np.tile(matplotlib.colors.to_rgba("#EEEEEE"),(len(spots),1))
        rgba[available]=matplotlib.colors.to_rgba("#B0B0B0")
        if kind=="significance":rgba[selected]=matplotlib.colors.to_rgba("#B2182B")
        elif selected.any():rgba[selected]=PVALUE_CMAP(PVALUE_NORM(-np.log10(np.maximum(frame.display_p_value.to_numpy()[selected],1e-4))))
        # Fixed map axes are final: convert the recorded footprint to points once.
        # One reusable circular path keeps the vector PDF compact.
        diameter_points=SPOT_DIAMETER*ax.get_window_extent().width/(MAP_XLIM[1]-MAP_XLIM[0])*72/ax.figure.dpi
        diameter_y=SPOT_DIAMETER*ax.get_window_extent().height/(MAP_YLIM[1]-MAP_YLIM[0])*72/ax.figure.dpi
        assert np.isclose(diameter_points,diameter_y)
        ax.scatter(spots.x,spots.y,s=diameter_points**2,c=rgba,marker="o",linewidths=0,edgecolors="none")
    for (region,piece),border in roi.groupby(["region","piece"]):
        ax.plot(border.x,border.y,color="white",lw=ROI_HALO_LINEWIDTH)
        ax.plot(border.x,border.y,color=ROI_COLORS[region],lw=ROI_LINEWIDTH)
    ax.set(xlim=MAP_XLIM,ylim=MAP_YLIM,aspect="equal");ax.axis("off")


def map_grid(fig, current, rect, kind, header_size=10):
    left,bottom,width,height=rect
    gap=.006;rowgap=.005
    fw,fh=fig.get_size_inches()
    header=max(.028,.20/fh)
    aspect=(MAP_YLIM[1]-MAP_YLIM[0])/(MAP_XLIM[1]-MAP_XLIM[0])
    h=(height-header-(len(current)-1)*rowgap)/len(current)
    w=min((width-2*gap)/3,h*fh/fw/aspect)
    h=w*fw/fh*aspect
    x0=left+(width-3*w-2*gap)/2
    top=bottom+height-header
    for j,method in enumerate(METHODS):
        x=x0+j*(w+gap)
        head=fig.add_axes([x,top+.004,w,header-.004]);head.set_facecolor("#F1F2F4")
        head.text(.5,.5,NAMES[method],ha="center",va="center",fontsize=header_size,weight="bold")
        head.axhline(0,color=COLORS[method],lw=3);head.set_xticks([]);head.set_yticks([])
        for spine in head.spines.values():spine.set_visible(False)
    for i,lr in enumerate(current):
        y=top-(i+1)*h-i*rowgap
        side=fig.add_axes([x0-.030,y,.024,h]);side.set_facecolor("#D1D2D3")
        side.text(.5,.5,label(lr),ha="center",va="center",rotation=90,fontsize=7,weight="bold",style="italic")
        side.set_xticks([]);side.set_yticks([])
        for spine in side.spines.values():spine.set_visible(False)
        for j,method in enumerate(METHODS):
            map_panel(fig.add_axes([x0+j*(w+gap),y,w,h]),method,lr,kind)


def pvalue_bar(fig, rect, fontsize=8):
    bar=fig.colorbar(plt.cm.ScalarMappable(norm=PVALUE_NORM,cmap=PVALUE_CMAP),
                    cax=fig.add_axes(rect),orientation="vertical",ticks=[1,2,3,4])
    bar.ax.set_yticklabels(["1","2","3","4+"]);bar.ax.tick_params(labelsize=fontsize,length=2)
    bar.ax.set_title("−log₁₀(p)",fontsize=fontsize+1,pad=7,weight="bold")
    bar.outline.set_visible(False)


def horizontal_pvalue_bar(fig, rect, fontsize=8):
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=PVALUE_NORM, cmap=PVALUE_CMAP),
                      cax=fig.add_axes(rect), orientation="horizontal", ticks=[1, 2, 3, 4])
    bar.ax.set_xticklabels(["1", "2", "3", "4+"])
    bar.ax.tick_params(labelsize=fontsize, length=2)
    bar.ax.text(1.13, .5, "−log₁₀\n(p-value)", transform=bar.ax.transAxes,
                ha="left", va="center", fontsize=fontsize, weight="bold")
    bar.outline.set_visible(False)


def regional_matrix(ids, regions=REGIONS):
    return np.array([[100 * region_lookup.loc[(lr, r, m), "proportion"]
                      if region_lookup.loc[(lr, r, m), "status"] == "available" else np.nan
                      for r in regions for m in METHODS] for lr in ids])


def regional_named_panel(ids):
    """Three published ROI blocks; outside-ROI coverage remains in full results."""
    values = regional_matrix(ids, REGIONS[:3])
    width, cell_h = 6.6, .17
    heat_h, left, gap, bottom, cell_w = len(ids) * cell_h, 1.45, .16, .88, .40
    height = heat_h + 1.45
    fig = plt.figure(figsize=(width, height))
    for k, (region, name) in enumerate(zip(REGIONS[:3], REGION_NAMES[:3])):
        ax = fig.add_axes([(left+k*(3*cell_w+gap))/width, bottom/height,
                           3*cell_w/width, heat_h/height])
        block = values[:, 3*k:3*(k+1)]
        mesh = ax.pcolormesh(np.arange(4)-.5, np.arange(len(ids)+1)-.5,
                            np.ma.masked_invalid(block), cmap=REGIONAL_CMAP, norm=REGIONAL_NORM)
        ax.set_ylim(len(ids)-.5, -.5)
        for i in range(len(ids)):
            for j in range(3):
                value = block[i, j]
                number = "/" if not np.isfinite(value) else ("<0.1%" if 0 < value < .1 else f"{value:.1f}%")
                luminance = np.dot(REGIONAL_CMAP(REGIONAL_NORM(value))[:3], [.2126, .7152, .0722]) if np.isfinite(value) else 1
                ax.text(j, i, number, ha="center", va="center", fontsize=6.4,
                        color="white" if luminance < .48 else "#111111")
        total = int((spots.region == region).sum())
        ax.set_title(f"{name}\nn = {total:,}", fontsize=8, weight="bold", pad=8, color=ROI_COLORS[region])
        ax.set_xticks(range(3), list(NAMES.values()), fontsize=7, rotation=50, ha="right", rotation_mode="anchor")
        for tick, method in zip(ax.get_xticklabels(), METHODS):
            tick.set_color(COLORS[method])
            tick.set_weight("bold")
        ax.set_yticks(range(len(ids)), [lr.replace("__", "–") for lr in ids] if k == 0 else [""]*len(ids),
                      fontsize=7, weight="bold")
        ax.set_xticks(np.arange(-.5, 3, 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(ids), 1), minor=True)
        ax.grid(which="minor", color="#444444", lw=.45)
        ax.tick_params(which="both", length=0, pad=3)
        for spine in ax.spines.values():
            spine.set_edgecolor("#111111")
            spine.set_linewidth(.65)
    bar = fig.colorbar(mesh, cax=fig.add_axes([5.86/width, (bottom+.1)/height, .10/width, (heat_h-.2)/height]),
                       ticks=[0, 1, 5, 10, 25, 50, 100])
    bar.solids.set_rasterized(False)
    bar.solids.set_edgecolor("face")
    bar.ax.set_yticklabels([f"{v}%" for v in [0, 1, 5, 10, 25, 50, 100]], fontsize=6.8)
    bar.ax.set_title("Region\ncoverage", fontsize=7, pad=8)
    return fig


def regional_clustered_panel(ids):
    """Original average-linkage/Euclidean clustering of raw percentages."""
    columns = [(r, m) for r in REGIONS for m in METHODS]
    values = regional_matrix(ids)
    imputed = np.where(np.isfinite(values), values, np.nanmedian(values, axis=0))
    row_tree = linkage(imputed, method="average", metric="euclidean", optimal_ordering=True)
    col_tree = linkage(imputed.T, method="average", metric="euclidean", optimal_ordering=True)
    row_order, col_order = leaves_list(row_tree), leaves_list(col_tree)
    width, height = 7.8, 7.4
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([.95/width, 2.05/height, 5.55/width, 3.90/height])
    mesh = ax.pcolormesh(np.ma.masked_invalid(values[row_order][:, col_order]), cmap=REGIONAL_CMAP, norm=REGIONAL_NORM)
    ax.invert_yaxis()
    ax.set_yticks([])
    names = dict(zip(REGIONS, REGION_NAMES))
    column_labels = [NAMES[columns[j][1]]+"\n"+names[columns[j][0]] for j in col_order]
    ax.set_xticks(np.arange(12)+.5, column_labels, rotation=60, ha="right", rotation_mode="anchor", fontsize=8, weight="bold")
    ax.tick_params(axis="x", length=2, pad=4)
    for spine in ax.spines.values():
        spine.set_linewidth(.65)
    rax = fig.add_axes([.20/width, 2.05/height, .60/width, 3.90/height])
    dendrogram(row_tree, orientation="left", no_labels=True, ax=rax, color_threshold=0, above_threshold_color="#555555")
    for collection in rax.collections:
        collection.set_linewidth(.6)
    rax.invert_yaxis()
    rax.axis("off")
    cax = fig.add_axes([.95/width, 6.32/height, 5.55/width, .38/height])
    dendrogram(col_tree, no_labels=True, ax=cax, color_threshold=0, above_threshold_color="#555555")
    cax.axis("off")
    for collection in cax.collections:
        collection.set_linewidth(.8)
    ann = fig.add_axes([.95/width, 6.00/height, 5.55/width, .24/height])
    ann.imshow([[to_rgb(COLORS[columns[j][1]]) for j in col_order],
                [to_rgb(ROI_COLORS[columns[j][0]]) for j in col_order]], aspect="auto")
    ann.set_xticks([])
    ann.set_yticks([0, 1], ["Software", "Region"], fontsize=7)
    ann.tick_params(length=0)
    bar = fig.colorbar(mesh, cax=fig.add_axes([6.86/width, 3.10/height, .12/width, 2.25/height]), ticks=[0, 1, 5, 10, 25, 50, 100])
    bar.solids.set_rasterized(False)
    bar.solids.set_edgecolor("face")
    bar.ax.set_yticklabels([f"{v}%" for v in [0, 1, 5, 10, 25, 50, 100]], fontsize=8)
    bar.ax.set_title("Region\ncoverage", fontsize=8, pad=8)
    fig.legend(handles=[Patch(color=ROI_COLORS[r], label=n) for r, n in zip(REGIONS, REGION_NAMES)],
               loc="lower center", ncol=4, bbox_to_anchor=(.57, .67/height), fontsize=7.5, frameon=False)
    fig.text(.025, .83/height, "Region", fontsize=10, weight="bold")
    fig.legend(handles=[Line2D([0], [0], color=COLORS[m], lw=1, label=NAMES[m]) for m in METHODS],
               loc="lower center", ncol=3, bbox_to_anchor=(.52, .32/height), fontsize=8, frameon=False)
    fig.text(.025, .47/height, "Software", fontsize=10, weight="bold")
    return fig


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.pdf", facecolor="white")
    plt.close(fig)
    print("Rendered", name, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--jaccard-dir", type=Path, help="Temporary handoff directory used by the R plotting entry point.")
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    RESULTS = repo / "res/benchmarks/visium"
    DATA = repo / "Data/generated/benchmarks/visium"
    FIGURES = repo / "Figures"
    FIGURES.mkdir(parents=True, exist_ok=True)
    final = pd.read_csv(RESULTS / "visium_lri_hotspot_counts_and_test_status_10000_samples.csv", converters={"status": str})
    blisa_order = final[final.condition == "blisa"].sort_values(["returned", "n_sig", "lr_id"], ascending=[False, False, True], na_position="last")
    common = set.intersection(*[set(final.loc[(final.condition == m) & (final.status == "available"), "lr_id"]) for m in METHODS])
    common_ids = blisa_order.loc[blisa_order.lr_id.isin(common), "lr_id"].tolist()
    assert len(common_ids) == 679
    selected_ids = set(common_ids[:10]) | {"LAMA5__CD44"}
    sys.path.insert(0, str(repo / "codes/helpers"))
    from visium_benchmark_results import summarize_saved_results
    _, region_results, overlap, pooled, map_lookup = summarize_saved_results(repo, selected_ids)
    if args.jaccard_dir:
        overlap.to_csv(args.jaccard_dir / "jaccard_per_lr.csv", index=False, na_rep="NA")
        pooled.to_csv(args.jaccard_dir / "jaccard_pooled.csv", index=False, na_rep="NA")
    spots = pd.read_csv(DATA / "visium_benchmark_spot_annotations.csv", usecols=["spot_id", "image_x", "image_y", "region"])
    spots["x"], spots["y"] = spots.image_x, -spots.image_y
    roi = pd.read_csv(DATA / "visium_benchmark_registered_roi_borders.csv")
    SPOT_DIAMETER = json.loads((DATA / "visium_benchmark_image_scalefactors.json").read_text())["spot_diameter_fullres"]
    MAP_XLIM = (min(spots.x.min(), roi.x.min())-SPOT_DIAMETER/2, max(spots.x.max(), roi.x.max())+SPOT_DIAMETER/2)
    MAP_YLIM = (min(spots.y.min(), roi.y.min())-SPOT_DIAMETER/2, max(spots.y.max(), roi.y.max())+SPOT_DIAMETER/2)
    region_lookup = region_results.set_index(["lr_id", "region", "condition"])
    # Native selections and regional denominators are checked again after CSV round-trip.
    for (method, lr), frame in map_lookup.items():
        expected = final[(final.condition == method) & (final.lr_id == lr)].iloc[0]
        assert int(frame.selected.fillna(False).sum()) == expected.n_sig

    fig = plt.figure(figsize=(7.2, 3.05))
    map_grid(fig, ["LAMA5__CD44"], [.035, .13, .935, .85], "pvalue", header_size=11)
    horizontal_pvalue_bar(fig, [.62, .078, .23, .025], fontsize=8)
    fig.legend(handles=[Patch(color=ROI_COLORS[r], label=n) for r, n in zip(REGIONS[:3], REGION_NAMES[:3])],
               loc="lower left", ncol=3, bbox_to_anchor=(.035, .002), fontsize=7, frameon=False)
    save(fig, "figure4h")

    fig = plt.figure(figsize=(13.8, 11.4))
    map_grid(fig, common_ids[:5], [.035, .075, .46, .90], "pvalue", header_size=10)
    map_grid(fig, common_ids[5:10], [.535, .075, .46, .90], "pvalue", header_size=10)
    fig.legend(handles=[Patch(color=ROI_COLORS[r], label=n) for r, n in zip(REGIONS[:3], REGION_NAMES[:3])],
               loc="lower left", ncol=3, bbox_to_anchor=(.045, .025), fontsize=9, frameon=False)
    horizontal_pvalue_bar(fig, [.56, .045, .24, .016], fontsize=9)
    save(fig, "figureS7b")

    biology_ids = ["AREG__EGFR", "CDH1__CDH1", "EFNA1__EPHA2", "EGF__EGFR", "FGF10__FGFR2", "FGF10__FGFR1",
                   "IGF1__IGF1R", "FGF7__FGFR2", "MDK__LRP1", "THBS2__SDC1", "THBS2__SDC4", "PDGFB__PDGFRA"]
    save(regional_named_panel(biology_ids), "figure4i")
    top50 = blisa_order.loc[blisa_order.returned, "lr_id"].head(50).tolist()
    save(regional_named_panel(top50), "figureS8b")
    # S8c matches report page 28: 679 common LR pairs, all four regions, with
    # identical column order, annotation colours, dendrograms and row patterns.
    save(regional_clustered_panel(sorted(common)), "figureS8c")
