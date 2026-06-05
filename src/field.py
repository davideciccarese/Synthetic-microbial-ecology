"""
field.py
========

A scalar diffusion field on a fine regular grid.

Each chemical the cells read or write (primary substrate, a leaked metabolite,
an antibiotic) is one ScalarField. The cells sample the field under their body
to set growth, and add or remove mass at the same nodes.

Two boundary regimes cover the biology in the paper:
  reservoir (Dirichlet): the substrate is supplied from outside the colony, as
    on an agar plate, so the border is held at a fixed concentration.
  zero flux (Neumann): a secreted molecule stays in the domain, so it builds up
    locally near the producers. This is exactly the localising effect of spatial
    structure that the review stresses: in a structured habitat diffusion keeps
    secreted molecules concentrated next to the cells that make them.
"""

import numpy as np


class ScalarField:
    def __init__(self, N, dx, D, dt, c0=0.0, boundary='reservoir',
                 reservoir=1.0, decay=0.0, max_sub=0.22):
        self.N = N
        self.dx = dx
        self.D = D
        self.dt = dt
        self.boundary = boundary
        self.reservoir = reservoir
        self.decay = decay
        self.c = np.full((N, N), c0, dtype=np.float64)
        # explicit Euler stability: D * dt_sub / dx^2 < 0.25
        lim = max_sub * dx * dx / max(D, 1e-9)
        self.nsub = max(1, int(np.ceil(dt / lim)))
        self.dt_sub = dt / self.nsub

    def diffuse(self):
        c = self.c
        D, dx2 = self.D, self.dx * self.dx
        for _ in range(self.nsub):
            lap = np.zeros_like(c)
            lap[1:-1, 1:-1] = (c[2:, 1:-1] + c[:-2, 1:-1]
                               + c[1:-1, 2:] + c[1:-1, :-2]
                               - 4 * c[1:-1, 1:-1]) / dx2
            c += self.dt_sub * (D * lap - self.decay * c)
            if self.boundary == 'reservoir':
                c[0, :] = self.reservoir
                c[-1, :] = self.reservoir
                c[:, 0] = self.reservoir
                c[:, -1] = self.reservoir
            else:  # zero flux (Neumann): copy the neighbour row
                c[0, :] = c[1, :]
                c[-1, :] = c[-2, :]
                c[:, 0] = c[:, 1]
                c[:, -1] = c[:, -2]
        np.clip(c, 0.0, None, out=c)

    # -- node indices for a set of world points -----------------------
    def nodes(self, pts):
        """Map world points (..., 2) to clamped integer grid indices."""
        ix = np.clip(np.round(pts[..., 0] / self.dx).astype(int), 0, self.N - 1)
        iy = np.clip(np.round(pts[..., 1] / self.dx).astype(int), 0, self.N - 1)
        return ix, iy

    def sample_cells(self, spine_pts):
        """Mean field value under each cell, shape (n,). spine_pts is (n,k,2)."""
        ix, iy = self.nodes(spine_pts)
        return self.c[ix, iy].mean(axis=1)

    def add_cells(self, spine_pts, per_cell):
        """Add per_cell (n,) mass to the nodes under each cell, split over k."""
        ix, iy = self.nodes(spine_pts)
        k = spine_pts.shape[1]
        share = (per_cell / k)[:, None] * np.ones((1, k))
        np.add.at(self.c, (ix.ravel(), iy.ravel()), share.ravel())
        np.clip(self.c, 0.0, None, out=self.c)
