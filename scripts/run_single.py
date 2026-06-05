"""
run_single.py
=============

Run one interaction and write its own GIF, in strain or lineage colours.

Usage:
    python scripts/run_single.py competition
    python scripts/run_single.py mutualism_active --mode lineage

Interaction names are the short row tags:
    commensalism, public_good, mutualism_passive, mutualism_active,
    competition, amensalism, predation
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from config import Config
import interactions as I
import sim
import render

BY_ROW = {x.row: x for x in I.ALL}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("interaction", choices=list(BY_ROW))
    ap.add_argument("--mode", default="strain", choices=["strain", "lineage"])
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    cfg = Config()
    inter = BY_ROW[args.interaction]
    frames, _ = sim.run(inter, cfg)

    here = os.path.dirname(__file__)
    figdir = os.path.abspath(os.path.join(here, "..", "figures"))
    os.makedirs(figdir, exist_ok=True)
    out = args.out or os.path.join(figdir, f"{args.interaction}_{args.mode}.gif")
    render.single_gif(inter, frames, cfg, out, mode=args.mode)
    print("wrote", out, "with", frames[-1].seg.shape[0], "cells")


if __name__ == "__main__":
    main()
