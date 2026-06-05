"""
render.py
=========

Turn snapshots into figures and GIFs.

Cells are drawn as capsule contours: the true spherocylinder outline of each
cell (a stadium shape) as a filled polygon with a thin dark edge, so individual
cells read crisply and, when left unfilled, an underlying field shows through.
Two colour modes for the cell panels:
  strain   blue for strain A, green for strain B
  lineage  a fixed colour per founding lineage, so sectors read as solid stripes

Three multipanel builders lay the seven interactions on one grid:
  multipanel_gif          cells coloured by strain or lineage
  multipanel_fields_gif   the representative diffusion field each colony shapes,
                          as a heatmap with the cells drawn as contours on top
"""

import numpy as np
import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

# substrate colour family: dark (low) -> bright green (high), so the "bright =
# high" reading is consistent with magma (metabolite) and hot (antibiotic)
_GREEN_GLOW = LinearSegmentedColormap.from_list(
    "green_glow", ["#021a0e", "#0a5a30", "#28c172", "#c9ffd9"])
try:
    matplotlib.colormaps.register(_GREEN_GLOW)
except Exception:
    pass


A_COLOR = "#2563eb"   # strain A, blue
B_COLOR = "#16a34a"   # strain B, green
EDGE = "#0f172a"      # dark cell outline for the capsule contour


def capsule_polys(seg, R, ncap=8):
    """Build a closed spherocylinder outline (capsule) for every cell.

    seg is (n, 2, 2), the two spine endpoints p and q. Each cell becomes a
    polygon: a rounded cap of ncap points at each end joined by the two straight
    flanks. Returned as an (n, 2*ncap, 2) array, ready for a PolyCollection.
    Drawing the true cell contour (rather than a thick line) gives crisp cell
    boundaries and, when unfilled, lets an underlying field show through.
    """
    if seg.shape[0] == 0:
        return np.zeros((0, 2 * ncap, 2))
    p = seg[:, 0, :]
    q = seg[:, 1, :]
    d = q - p
    phi = np.arctan2(d[:, 1], d[:, 0])
    aq = phi[:, None] + np.linspace(np.pi / 2, -np.pi / 2, ncap)[None, :]
    ap = phi[:, None] + np.linspace(-np.pi / 2, -3 * np.pi / 2, ncap)[None, :]
    arc_q = np.stack([q[:, 0, None] + R * np.cos(aq),
                      q[:, 1, None] + R * np.sin(aq)], axis=2)
    arc_p = np.stack([p[:, 0, None] + R * np.cos(ap),
                      p[:, 1, None] + R * np.sin(ap)], axis=2)
    return np.concatenate([arc_q, arc_p], axis=1)


def lineage_colors(max_lin, seed=0):
    rng = np.random.default_rng(seed)
    base = plt.get_cmap("twilight")(np.linspace(0, 1, max_lin + 1))
    rng.shuffle(base)
    return base


def _linewidth_points(ax, R, fig):
    """Convert a rod half width in data units to a line width in points."""
    # data-to-display scale on x, then display points
    bbox = ax.get_window_extent()
    xlim = ax.get_xlim()
    px_per_data = bbox.width / (xlim[1] - xlim[0])
    pts_per_px = 72.0 / fig.dpi
    return 2 * R * px_per_data * pts_per_px


def draw_panel(ax, snap, view, mode="strain", lin_colors=None, fig=None):
    ax.clear()
    ax.set_xlim(view[0], view[1])
    ax.set_ylim(view[0], view[1])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#d4d4d8")
        s.set_linewidth(0.8)
    if snap.seg.shape[0] == 0:
        return
    if mode == "strain":
        cols = np.where(snap.sp == 0, A_COLOR, B_COLOR)
    else:
        cols = lin_colors[snap.lin % len(lin_colors)]
    lw = _linewidth_points(ax, snap.R, fig) if fig is not None else 2.0
    lc = LineCollection(snap.seg, colors=cols, linewidths=lw,
                        capstyle="round", joinstyle="round")
    ax.add_collection(lc)


def _fig_to_image(fig):
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    return Image.fromarray(buf[..., :3].copy())


def _panel_colors(snap, mode, lin_colors):
    if mode == "strain":
        return np.where(snap.sp == 0, A_COLOR, B_COLOR)
    return lin_colors[snap.lin % len(lin_colors)]


def _setup_panel(ax, view, title, caption):
    """Style a panel once: limits, frame, title and caption do not change."""
    ax.set_xlim(view[0], view[1])
    ax.set_ylim(view[0], view[1])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#d4d4d8")
        s.set_linewidth(0.9)
    ax.set_title(textwrap.fill(title, width=24), fontsize=15,
                 color="#18181b", pad=7)
    ax.text(0.5, -0.05, caption, transform=ax.transAxes, ha="center",
            va="top", fontsize=13, color="#52525b")


def multipanel_gif(results, cfg, path, mode="strain", fps=8, title=None):
    """results: list of (interaction, frames). Writes an animated GIF.

    Layout is built once and never recomputed, so every panel stays locked in
    place. Only the rod geometry inside each panel changes between frames.
    """
    ncols = 4
    nrows = 2
    n_frames = cfg.n_frames

    # one common, fixed view box big enough for the largest final colony
    c = cfg.center
    maxr = 0.0
    for _, frames in results:
        seg = frames[-1].seg
        if seg.shape[0]:
            r = np.abs(seg.reshape(-1, 2) - c).max()
            maxr = max(maxr, r)
    half = maxr * 1.10 + 3
    view = (c - half, c + half)

    lin_colors = lineage_colors(cfg.n_seed + 5, seed=1) if mode == "lineage" else None

    fig, axes = plt.subplots(nrows, ncols, figsize=(17.0, 9.4), dpi=120)
    fig.patch.set_facecolor("white")
    axes = axes.ravel()
    legend_ax = axes[-1]

    # fixed margins, set once. No tight_layout inside the loop, which is what
    # was nudging the panels frame to frame.
    fig.subplots_adjust(left=0.012, right=0.988, top=0.875, bottom=0.05,
                        wspace=0.10, hspace=0.32)
    if title:
        fig.suptitle(title, fontsize=21, color="#18181b", y=0.975)

    # style each interaction panel once and attach an empty, persistent
    # PolyCollection of capsule contours we update in place every frame
    collections = []
    for k, (inter, _frames) in enumerate(results):
        ax = axes[k]
        caption = f"effect on (A, B): ({inter.signs[0]}, {inter.signs[1]})"
        _setup_panel(ax, view, inter.name, caption)
        pc = PolyCollection([], edgecolors=EDGE, linewidths=0.3)
        ax.add_collection(pc)
        collections.append(pc)
    _draw_legend(legend_ax, cfg, mode)

    fig.canvas.draw()

    imgs = []
    for fi in range(n_frames):
        for k, (inter, frames) in enumerate(results):
            snap = frames[fi]
            pc = collections[k]
            pc.set_verts(capsule_polys(snap.seg, snap.R))
            pc.set_facecolor(_panel_colors(snap, mode, lin_colors))
        imgs.append(_fig_to_image(fig))

    plt.close(fig)
    dur = int(1000 / fps)
    imgs[0].save(path, save_all=True, append_images=imgs[1:],
                 duration=dur, loop=0, optimize=True)
    return path


def multipanel_fields_gif(results, cfg, path, fps=8, title=None):
    """Seven-panel GIF of the field each colony shapes, cells drawn as contours.

    results: list of (interaction, frames, field_label, field_arrays, cmap), one
    per interaction. field_arrays is a list of N x N grids (one per frame) for
    the representative molecule, and cmap sets its colour family. Each panel
    shows that field as a heatmap (bright = high) with the colony outline on top.
    """
    ncols, nrows = 4, 2
    n_frames = cfg.n_frames
    c = cfg.center
    maxr = 0.0
    for _, frames, _, _, _ in results:
        seg = frames[-1].seg
        if seg.shape[0]:
            maxr = max(maxr, np.abs(seg.reshape(-1, 2) - c).max())
    half = maxr * 1.10 + 3
    view = (c - half, c + half)
    extent = [0, cfg.N * cfg.dx, 0, cfg.N * cfg.dx]

    fig, axes = plt.subplots(nrows, ncols, figsize=(17.0, 9.4), dpi=120)
    fig.patch.set_facecolor("white")
    axes = axes.ravel()
    legend_ax = axes[-1]
    fig.subplots_adjust(left=0.012, right=0.988, top=0.875, bottom=0.05,
                        wspace=0.10, hspace=0.32)
    if title:
        fig.suptitle(title, fontsize=21, color="#18181b", y=0.975)

    ims, pcs = [], []
    for k, (inter, frames, label, arrs, cmap) in enumerate(results):
        ax = axes[k]
        _setup_panel(ax, view, inter.name, label)
        vmax = max(max(a.max() for a in arrs), 1e-6)
        im = ax.imshow(arrs[0].T, origin="lower", extent=extent, cmap=cmap,
                       vmin=0.0, vmax=vmax, interpolation="bilinear", zorder=0)
        ax.set_xlim(view[0], view[1])
        ax.set_ylim(view[0], view[1])
        pc = PolyCollection([], facecolors="none", linewidths=0.6, zorder=2)
        ax.add_collection(pc)
        ims.append(im)
        pcs.append(pc)
    _draw_legend(legend_ax, cfg, "fields")
    fig.canvas.draw()

    imgs = []
    for fi in range(n_frames):
        for k, (inter, frames, label, arrs, cmap) in enumerate(results):
            snap = frames[fi]
            ims[k].set_data(arrs[fi].T)
            pcs[k].set_verts(capsule_polys(snap.seg, snap.R))
            pcs[k].set_edgecolor(np.where(snap.sp == 0, A_COLOR, B_COLOR))
        imgs.append(_fig_to_image(fig))
    plt.close(fig)
    imgs[0].save(path, save_all=True, append_images=imgs[1:],
                 duration=int(1000 / fps), loop=0, optimize=True)
    return path


def _draw_legend(ax, cfg, mode):
    ax.clear()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    if mode == "fields":
        ax.text(0.02, 0.98, "How to read the panels", fontsize=16.5,
                color="#18181b", va="top", weight="bold")
        ax.text(0.02, 0.87, "Cells (outline by strain)", fontsize=13.5,
                color="#3f3f46", va="center")
        for col, lab, yy in [("#2563eb", "Strain A", 0.79),
                             ("#16a34a", "Strain B", 0.72)]:
            ax.add_patch(matplotlib.patches.Ellipse(
                (0.09, yy), 0.10, 0.045, facecolor="none", edgecolor=col,
                lw=2.0, transform=ax.transAxes))
            ax.text(0.20, yy, lab, fontsize=12.5, color="#3f3f46", va="center")
        ax.text(0.02, 0.64, "Heatmap = molecule concentration\n(dark = low, bright = high):",
                fontsize=11.5, color="#3f3f46", va="top", linespacing=1.35)
        families = [("green_glow", "Substrate (resource)", 0.50),
                    ("magma", "Metabolite / public good", 0.40),
                    ("hot", "Antibiotic (toxin)", 0.30)]
        grad = np.linspace(0, 1, 256).reshape(1, -1)
        for cmap, lab, yy in families:
            cax = ax.inset_axes([0.04, yy, 0.34, 0.032])
            cax.imshow(grad, aspect="auto", cmap=cmap)
            cax.set_xticks([])
            cax.set_yticks([])
            ax.text(0.42, yy + 0.016, lab, fontsize=11, color="#52525b",
                    va="center")
        notes = (
            "Substrate is supplied from the border and\n"
            "drawn down under the colony. Secreted\n"
            "molecules stay local, building up next to\n"
            "their producers. Each panel scaled to its\n"
            "own maximum."
        )
        ax.text(0.02, 0.21, notes, fontsize=10.5, color="#52525b", va="top",
                linespacing=1.4)
        return

    ax.text(0.02, 0.98, "How to read the panels", fontsize=16.5,
            color="#18181b", va="top", weight="bold")

    if mode == "lineage":
        # a colour bar standing in for the founder lineage palette: each founder
        # cell keeps its own hue, so a clonal sector reads as one solid colour
        ax.text(0.02, 0.88, "Colour = founder lineage", fontsize=14,
                color="#3f3f46", va="center")
        cax = ax.inset_axes([0.04, 0.74, 0.66, 0.055])
        grad = np.linspace(0, 1, 256).reshape(1, -1)
        cax.imshow(grad, aspect="auto", cmap="twilight")
        cax.set_xticks([0, 255])
        cax.set_xticklabels(["founder 1", f"founder {cfg.n_seed}"], fontsize=10)
        cax.set_yticks([])
        cax.tick_params(length=0)
        y = 0.62
    else:
        lines = [("#2563eb", "Strain A"), ("#16a34a", "Strain B")]
        y = 0.88
        for col, lab in lines:
            ax.add_line(plt.Line2D([0.04, 0.16], [y, y], color=col, lw=6,
                                   solid_capstyle="round"))
            ax.text(0.20, y, lab, fontsize=14, color="#3f3f46", va="center")
            y -= 0.095
        y -= 0.01

    notes = (
        "All panels share one founder disk, identical\n"
        "rod mechanics and substrate supply. Only the\n"
        "chemical coupling differs, so every spatial\n"
        "outcome is emergent.\n\n"
        "Sign key (effect of partner on focal strain):\n"
        "  +  promoted    -  reduced    0  no effect\n\n"
        "Substrate S is supplied from the border, like\n"
        "an agar plate. Secreted molecules stay local,\n"
        "building up by their producers. Only the rim\n"
        "grows; the buried core stays frozen."
    )
    ax.text(0.02, y, notes, fontsize=11, color="#52525b", va="top",
            linespacing=1.4)


def _overlay_segments(ax, view):
    """Thin strain-coloured colony overlay drawn on top of a field heatmap."""
    ax.set_xlim(view[0], view[1])
    ax.set_ylim(view[0], view[1])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    lc = LineCollection([], capstyle="round", joinstyle="round",
                        linewidths=1.4, alpha=0.9)
    ax.add_collection(lc)
    return lc


def fields_gif(inter, frames, field_hist, cfg, path, fps=8,
               field_specs=None):
    """Animate the colony next to the diffusion fields it shapes.

    Left panel: the cells, coloured by strain. Each further panel is one
    chemical field (substrate, secreted metabolite, antibiotic) as a heatmap
    with the colony outline on top, so the localised depletion and build up of
    each molecule under the colony is visible. This is the spatial structure
    the review emphasises: secreted molecules stay next to their producers.
    """
    if field_specs is None:
        field_specs = [(k, f"Field {k}", "magma") for k in field_hist]
    keys = [k for k, _, _ in field_specs if k in field_hist]
    n_fields = len(keys)
    ncols = 1 + n_fields

    c = cfg.center
    seg = frames[-1].seg
    maxr = np.abs(seg.reshape(-1, 2) - c).max() if seg.shape[0] else 30
    half = maxr * 1.12 + 3
    view = (c - half, c + half)
    extent = [0, cfg.N * cfg.dx, 0, cfg.N * cfg.dx]

    # fixed colour scale per field across all frames
    vmax = {}
    for k in keys:
        m = max(arr.max() for arr in field_hist[k])
        vmax[k] = max(m, 1e-6)

    fig, axes = plt.subplots(1, ncols, figsize=(4.7 * ncols, 5.6), dpi=120)
    if ncols == 1:
        axes = [axes]
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.04,
                        wspace=0.12)
    fig.suptitle(f"{inter.name}: cells and the fields they shape",
                 fontsize=16, color="#18181b", y=0.97)

    ax_cells = axes[0]
    ax_cells.set_title("Cells (blue A, green B)", fontsize=13.5,
                       color="#18181b", pad=7)
    for s in ax_cells.spines.values():
        s.set_color("#d4d4d8")
    ax_cells.set_xlim(view[0], view[1])
    ax_cells.set_ylim(view[0], view[1])
    ax_cells.set_aspect("equal")
    ax_cells.set_xticks([])
    ax_cells.set_yticks([])
    cells_pc = PolyCollection([], edgecolors=EDGE, linewidths=0.3)
    ax_cells.add_collection(cells_pc)

    titles = {k: t for k, t, _ in field_specs}
    cmaps = {k: cm for k, _, cm in field_specs}
    ims = {}
    overlays = {}
    for col_i, k in enumerate(keys, start=1):
        ax = axes[col_i]
        ax.set_title(titles[k], fontsize=13.5, color="#18181b", pad=7)
        im = ax.imshow(field_hist[k][0].T, origin="lower", extent=extent,
                       cmap=cmaps[k], vmin=0.0, vmax=vmax[k],
                       interpolation="bilinear")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=10)
        ims[k] = im
        ax.set_xlim(view[0], view[1])
        ax.set_ylim(view[0], view[1])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        pc = PolyCollection([], facecolors="none", linewidths=0.6, zorder=2)
        ax.add_collection(pc)
        overlays[k] = pc

    fig.canvas.draw()

    imgs = []
    for fi in range(len(frames)):
        snap = frames[fi]
        cols = np.where(snap.sp == 0, A_COLOR, B_COLOR)
        polys = capsule_polys(snap.seg, snap.R)
        cells_pc.set_verts(polys)
        cells_pc.set_facecolor(cols)
        for k in keys:
            ims[k].set_data(field_hist[k][fi].T)
            overlays[k].set_verts(polys)
            overlays[k].set_edgecolor(cols)
        imgs.append(_fig_to_image(fig))
    plt.close(fig)
    imgs[0].save(path, save_all=True, append_images=imgs[1:],
                 duration=int(1000 / fps), loop=0, optimize=True)
    return path


def single_gif(inter, frames, cfg, path, mode="strain", fps=8):
    c = cfg.center
    seg = frames[-1].seg
    maxr = np.abs(seg.reshape(-1, 2) - c).max() if seg.shape[0] else 30
    half = maxr * 1.10 + 3
    view = (c - half, c + half)
    lin_colors = lineage_colors(cfg.n_seed + 5, seed=1) if mode == "lineage" else None

    fig, ax = plt.subplots(figsize=(6.2, 6.6), dpi=130)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.03, right=0.97, top=0.93, bottom=0.03)
    caption = f"effect on (A, B): ({inter.signs[0]}, {inter.signs[1]})"
    _setup_panel(ax, view, inter.name, caption)
    ax.title.set_fontsize(15)
    pc = PolyCollection([], edgecolors=EDGE, linewidths=0.3)
    ax.add_collection(pc)
    fig.canvas.draw()

    imgs = []
    for fr in frames:
        pc.set_verts(capsule_polys(fr.seg, fr.R))
        pc.set_facecolor(_panel_colors(fr, mode, lin_colors))
        imgs.append(_fig_to_image(fig))
    plt.close(fig)
    imgs[0].save(path, save_all=True, append_images=imgs[1:],
                 duration=int(1000 / fps), loop=0, optimize=True)
    return path
