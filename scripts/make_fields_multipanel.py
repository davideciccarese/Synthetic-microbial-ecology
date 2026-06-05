"""
make_fields_multipanel.py
=========================

Seven-panel companion to the strain and lineage multipanels. Each panel shows
the representative diffusion field the colony shapes (the substrate it depletes,
or the metabolite, public good or antibiotic it builds up) as a heatmap, with
the cells drawn as contours on top. Colour families code the molecule type:
substrate (resource) in green, metabolite or public good in magma, antibiotic in
hot.

Usage:
    python scripts/make_fields_multipanel.py [--out PATH] [--frames N]
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from config import Config
import interactions as I
import sim
import render


# colour family per molecule type
CMAP = {"S": "green_glow", "M": "magma", "Ma": "magma", "Mb": "magma",
        "P": "magma", "T": "hot"}

# representative field per interaction: (field key, short panel label)
REPRESENTATIVE = {
    "commensalism": ("M", "Metabolite M (leaked by A)"),
    "public_good": ("P", "Public good P (made by A)"),
    "mutualism_passive": ("Ma", "Metabolite from A (by-product)"),
    "mutualism_active": ("Ma", "Metabolite from A (obligate)"),
    "competition": ("S", "Substrate S (the only resource)"),
    "amensalism": ("T", "Antibiotic T (made by A)"),
    "predation": ("S", "Substrate S (prey resource)"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--frames", type=int, default=None)
    args = ap.parse_args()

    cfg = Config()
    if args.frames:
        cfg.n_frames = args.frames

    here = os.path.dirname(__file__)
    figdir = os.path.abspath(os.path.join(here, "..", "figures"))
    os.makedirs(figdir, exist_ok=True)
    out = args.out or os.path.join(figdir, "fields_multipanel.gif")

    results = []
    for inter in I.ALL:
        key, label = REPRESENTATIVE[inter.row]
        t0 = time.time()
        frames, _F, field_hist = sim.run_with_fields(inter, cfg, [key])
        results.append((inter, frames, label, field_hist[key], CMAP[key]))
        print(f"  {inter.row:20s} field={key:2s} "
              f"{frames[-1].seg.shape[0]:5d} cells {time.time() - t0:5.1f}s")

    render.multipanel_fields_gif(
        results, cfg, out,
        title="The fields each colony shapes (after Dolinsek et al. 2016)")
    print("wrote", out)


if __name__ == "__main__":
    main()
