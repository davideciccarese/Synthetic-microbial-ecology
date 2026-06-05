"""
make_multipanel.py
==================

Run all seven interaction archetypes from the same founder disk and write the
multipanel GIF. This is the hero figure of the repository.

Usage:
    python scripts/make_multipanel.py [--mode strain|lineage] [--out PATH]
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="strain", choices=["strain", "lineage"])
    ap.add_argument("--out", default=None)
    ap.add_argument("--frames", type=int, default=None)
    args = ap.parse_args()

    cfg = Config()
    if args.frames:
        cfg.n_frames = args.frames

    here = os.path.dirname(__file__)
    figdir = os.path.abspath(os.path.join(here, "..", "figures"))
    os.makedirs(figdir, exist_ok=True)
    out = args.out or os.path.join(figdir, f"interactions_{args.mode}.gif")

    results = []
    for inter in I.ALL:
        t0 = time.time()
        frames, _ = sim.run(inter, cfg)
        results.append((inter, frames))
        print(f"  {inter.row:22s} {frames[-1].seg.shape[0]:5d} cells "
              f"{time.time() - t0:5.1f}s")

    render.multipanel_gif(
        results, cfg, out, mode=args.mode,
        title="Microbial interaction archetypes in a range expanding rod colony "
              "(after Dolinsek et al. 2016, Fig. 1)")
    print("wrote", out)


if __name__ == "__main__":
    main()
