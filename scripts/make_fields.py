"""
make_fields.py
==============

Render a single interaction as the colony next to the diffusion fields it
shapes: the substrate it depletes, the metabolite or public good it builds up,
the antibiotic it spreads. Colour families code the molecule type (substrate in
green, metabolite or public good in magma, antibiotic in hot).

Usage:
    python scripts/make_fields.py [--row commensalism] [--out PATH] [--frames N]

Rows: commensalism, public_good, mutualism_passive, mutualism_active,
      competition, amensalism, predation
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


# (field key, label, colormap) per interaction
FIELD_SPECS = {
    "commensalism": [("S", "Substrate S (consumed by A)", "green_glow"),
                     ("M", "Metabolite M (leaked by A, eaten by B)", "magma")],
    "public_good": [("S", "Substrate S (shared)", "green_glow"),
                    ("P", "Public good P (made by A, used by B)", "magma")],
    "mutualism_passive": [("Ma", "By-product from A (feeds B)", "magma"),
                          ("Mb", "By-product from B (feeds A)", "magma")],
    "mutualism_active": [("Ma", "Metabolite from A (B needs it)", "magma"),
                         ("Mb", "Metabolite from B (A needs it)", "magma")],
    "competition": [("S", "Substrate S (the only resource)", "green_glow")],
    "predation": [("S", "Substrate S (prey resource)", "green_glow")],
    "amensalism": [("S", "Substrate S", "green_glow"),
                   ("T", "Antibiotic T (made by A, inhibits B)", "hot")],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--row", default="commensalism",
                    choices=list(FIELD_SPECS.keys()))
    ap.add_argument("--out", default=None)
    ap.add_argument("--frames", type=int, default=None)
    args = ap.parse_args()

    cfg = Config()
    if args.frames:
        cfg.n_frames = args.frames

    by_row = {inter.row: inter for inter in I.ALL}
    inter = by_row[args.row]
    specs = FIELD_SPECS[args.row]
    keys = [k for k, _, _ in specs]

    here = os.path.dirname(__file__)
    figdir = os.path.abspath(os.path.join(here, "..", "figures"))
    os.makedirs(figdir, exist_ok=True)
    out = args.out or os.path.join(figdir, f"fields_{args.row}.gif")

    t0 = time.time()
    frames, _F, field_hist = sim.run_with_fields(inter, cfg, keys)
    print(f"  {args.row:20s} {frames[-1].seg.shape[0]:5d} cells "
          f"{time.time() - t0:5.1f}s")

    render.fields_gif(inter, frames, field_hist, cfg, out, field_specs=specs)
    print("wrote", out)


if __name__ == "__main__":
    main()
