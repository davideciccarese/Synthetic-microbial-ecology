"""
sim.py
======

Wire a RodColony to a set of diffusion fields and one Interaction, advance the
system, and capture a snapshot per frame. A snapshot is just the rod geometry
plus strain and lineage labels, which is everything the renderer needs.

Step order each macro step:
  1. interaction.step  -> per cell elongation rate mu, and field edits
  2. elongate by mu
  3. divide cells that grew past threshold
  4. relax contacts (overdamped pushing and turning)
  5. remove any cells killed this step (predation, antibiotic)
  6. diffuse every field
"""

import numpy as np

from rods import RodColony


class Snapshot:
    __slots__ = ("seg", "sp", "lin", "R")

    def __init__(self, seg, sp, lin, R):
        self.seg = seg    # (n, 2, 2) segment endpoints for LineCollection
        self.sp = sp      # (n,) strain label
        self.lin = lin    # (n,) lineage id
        self.R = R        # rod half width, for line width in data units


def _snapshot(col, cfg):
    p, q = col.endpoints()
    seg = np.stack([p, q], axis=1)
    return Snapshot(seg.copy(), col.sp.copy(), col.lin.copy(), cfg.R)


def run(interaction, cfg):
    """Run one interaction to completion and return a list of Snapshot."""
    frames, F, _ = _run_core(interaction, cfg, capture_fields=())
    return frames, F


def run_with_fields(interaction, cfg, field_keys):
    """Like run, but also return per-frame copies of the named fields.

    Returns (frames, F, field_hist) where field_hist[key] is a list of N x N
    grid copies, one per saved frame, for heatmap panels of nutrients and
    secreted metabolites.
    """
    return _run_core(interaction, cfg, capture_fields=tuple(field_keys))


def _run_core(interaction, cfg, capture_fields=()):
    rng = np.random.default_rng(cfg.seed)
    col = RodColony(cfg, rng)
    frac = interaction.seed_frac if interaction.seed_frac is not None else cfg.frac_strain1
    col.seed_disk(cfg.center, cfg.center, cfg.seed_radius, cfg.n_seed, frac)
    F = interaction.fields(cfg)

    # settle the initial packing before growth so founders do not overlap
    col.relax(iters=12)

    field_hist = {k: [] for k in capture_fields}

    def grab_fields():
        for k in capture_fields:
            field_hist[k].append(F[k].c.copy())

    frames = [_snapshot(col, cfg)]
    grab_fields()
    for _ in range(cfg.n_frames - 1):
        for _ in range(cfg.steps_per_frame):
            mu = interaction.step(col, F, cfg, cfg.dt)
            # range expansion: only the rim grows. Gate elongation by how close
            # each cell is to the front (computed on the current cells).
            factor, _ = col.front_factor()
            mu = mu * factor
            col._frozen = None
            col.elongate(mu, cfg.dt)
            col.divide()
            # freeze the jammed core for the contact solve (recomputed after
            # division so the mask length matches the new cell count)
            _, frozen = col.front_factor()
            col._frozen = frozen
            col.relax()
            if not col.alive.all():
                col.compact()
            for f in F.values():
                f.diffuse()
        frames.append(_snapshot(col, cfg))
        grab_fields()
    return frames, F, field_hist
