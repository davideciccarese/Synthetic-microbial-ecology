"""
rods.py
=======

Spherocylinder individual based rod colony.

Each cell is a rod with hemispherical caps described by:
  center (x, y), orientation angle theta, cylinder length L, fixed half width R.
The two spine endpoints are center +/- (L/2) * (cos theta, sin theta).

Cells do four things:
  - elongate at a rate set externally (by the ecological interaction module),
  - divide along the long axis once they grow past a threshold length,
  - push and rotate their neighbours through overdamped contact relaxation,
  - carry a permanent strain label and a permanent lineage id.

The mechanics here are deliberately the minimal core shared by published off
lattice engines (CellModeller, iDynoMiCS, BSim). Everything ecological lives in
separate modules so the same mechanics can play out every interaction type.

This engine reuses the conventions validated in earlier work on rod shaped
cross feeding colonies (Ciccarese et al., The ISME Journal 16:1453-1463, 2022).
"""

import numpy as np
from scipy.spatial import cKDTree


# ----------------------------------------------------------------------
# Geometry helper: closest distance between two 2D segments
# ----------------------------------------------------------------------
def _segment_closest_points(p1, q1, p2, q2):
    """Closest points between segment p1-q1 and segment p2-q2 in 2D.

    Returns (c1, c2) the closest point on each segment. Standard clamped
    parametric solution (Ericson, Real Time Collision Detection).
    """
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    a = d1 @ d1
    e = d2 @ d2
    f = d2 @ r

    if a <= 1e-12 and e <= 1e-12:
        return p1, p2
    if a <= 1e-12:
        s = 0.0
        t = np.clip(f / e, 0.0, 1.0)
    else:
        c = d1 @ r
        if e <= 1e-12:
            t = 0.0
            s = np.clip(-c / a, 0.0, 1.0)
        else:
            b = d1 @ d2
            denom = a * e - b * b
            s = np.clip((b * f - c * e) / denom, 0.0, 1.0) if denom > 1e-12 else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t = 0.0
                s = np.clip(-c / a, 0.0, 1.0)
            elif t > 1.0:
                t = 1.0
                s = np.clip((b - c) / a, 0.0, 1.0)
    c1 = p1 + d1 * s
    c2 = p2 + d2 * t
    return c1, c2


def _segment_closest_batch(p1, q1, p2, q2):
    """Vectorized closest points between two stacks of segments (M, 2).

    Same clamped parametric solution as the scalar version, evaluated for all M
    candidate pairs at once. This is what makes the contact step fast enough to
    push a colony of a couple of thousand rods.
    """
    eps = 1e-9
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    a = np.einsum("ij,ij->i", d1, d1)
    e = np.einsum("ij,ij->i", d2, d2)
    f = np.einsum("ij,ij->i", d2, r)
    c = np.einsum("ij,ij->i", d1, r)
    b = np.einsum("ij,ij->i", d1, d2)
    a = np.maximum(a, eps)
    e = np.maximum(e, eps)
    denom = a * e - b * b
    s = np.where(denom > eps, (b * f - c * e) / np.where(denom > eps, denom, 1.0), 0.0)
    s = np.clip(s, 0.0, 1.0)
    t = (b * s + f) / e
    # t < 0 branch
    lo = t < 0.0
    t = np.where(lo, 0.0, t)
    s = np.where(lo, np.clip(-c / a, 0.0, 1.0), s)
    # t > 1 branch
    hi = t > 1.0
    t = np.where(hi, 1.0, t)
    s = np.where(hi, np.clip((b - c) / a, 0.0, 1.0), s)
    c1 = p1 + d1 * s[:, None]
    c2 = p2 + d2 * t[:, None]
    return c1, c2


class RodColony:
    """A growing colony of spherocylinder cells.

    Attributes that matter to the rest of the code, all numpy arrays of length n:
      x, y     cell centre coordinates
      th       orientation angle
      L        current cylinder length (caps add R at each end)
      sp       strain label (0 or 1)
      lin      permanent lineage id (daughters inherit the founder id)
      alive    boolean mask (predation can kill cells)
    """

    def __init__(self, cfg, rng):
        self.cfg = cfg
        self.rng = rng
        self.x = np.zeros(0)
        self.y = np.zeros(0)
        self.th = np.zeros(0)
        self.L = np.zeros(0)
        self.sp = np.zeros(0, dtype=np.int8)
        self.lin = np.zeros(0, dtype=np.int32)
        self.alive = np.zeros(0, dtype=bool)
        self._next_lin = 0

    # -- seeding ------------------------------------------------------
    def seed_disk(self, cx, cy, radius, n_seed, frac_strain1):
        """Seed a central mixed disk of founder cells (range expansion start)."""
        rng = self.rng
        rr = radius * np.sqrt(rng.random(n_seed))
        aa = rng.random(n_seed) * 2 * np.pi
        self.x = cx + rr * np.cos(aa)
        self.y = cy + rr * np.sin(aa)
        self.th = rng.random(n_seed) * 2 * np.pi
        self.L = np.full(n_seed, self.cfg.L_birth)
        self.sp = (rng.random(n_seed) < frac_strain1).astype(np.int8)
        self.lin = np.arange(n_seed, dtype=np.int32)
        self.alive = np.ones(n_seed, dtype=bool)
        self._next_lin = n_seed

    # -- convenience --------------------------------------------------
    @property
    def n(self):
        return self.x.size

    def endpoints(self):
        """Return spine endpoint arrays p (n,2) and q (n,2)."""
        hx = 0.5 * self.L * np.cos(self.th)
        hy = 0.5 * self.L * np.sin(self.th)
        p = np.column_stack((self.x - hx, self.y - hy))
        q = np.column_stack((self.x + hx, self.y + hy))
        return p, q

    def spine_samples(self, k=3):
        """Sample k points evenly along each cell spine, shape (n, k, 2).

        Used to read and write the diffusion fields under the cell body.
        """
        p, q = self.endpoints()
        ts = np.linspace(0.0, 1.0, k)
        pts = p[:, None, :] * (1 - ts)[None, :, None] + q[:, None, :] * ts[None, :, None]
        return pts

    # -- growth and division -----------------------------------------
    def elongate(self, mu, dt):
        """Elongate every cell by mu * dt (mu is a per cell length rate)."""
        self.L = self.L + np.maximum(mu, 0.0) * dt

    def divide(self):
        """Split cells that exceed L_div into two daughters end to end."""
        cfg = self.cfg
        ready = (self.L >= cfg.L_div) & self.alive
        if not ready.any():
            return
        # cap total population so runtime stays bounded
        room = cfg.max_cells - self.n
        idx = np.where(ready)[0]
        if room <= 0:
            # at the cap, stop dividing but let cells stay at L_div
            self.L[idx] = cfg.L_div
            return
        if idx.size > room:
            idx = idx[:room]

        gap = 2 * cfg.R
        new_L = 0.5 * (self.L[idx] - gap)
        new_L = np.maximum(new_L, cfg.L_birth * 0.6)
        ux = np.cos(self.th[idx])
        uy = np.sin(self.th[idx])
        off = 0.5 * (new_L + gap)

        # daughter A keeps the slot, daughter B is appended
        ax = self.x[idx] - off * ux
        ay = self.y[idx] - off * uy
        bx = self.x[idx] + off * ux
        by = self.y[idx] + off * uy
        noise = self.cfg.div_angle_noise
        tha = self.th[idx] + self.rng.normal(0, noise, idx.size)
        thb = self.th[idx] + self.rng.normal(0, noise, idx.size)

        self.x[idx] = ax
        self.y[idx] = ay
        self.th[idx] = tha
        self.L[idx] = new_L

        self.x = np.concatenate([self.x, bx])
        self.y = np.concatenate([self.y, by])
        self.th = np.concatenate([self.th, thb])
        self.L = np.concatenate([self.L, new_L])
        self.sp = np.concatenate([self.sp, self.sp[idx]])
        self.lin = np.concatenate([self.lin, self.lin[idx]])
        self.alive = np.concatenate([self.alive, np.ones(idx.size, dtype=bool)])

    # -- range expansion front ---------------------------------------
    def front_factor(self):
        """Per cell growth multiplier in [0, 1] and a deep-core freeze mask.

        A cell at the colony edge has neighbours only on one side, so the unit
        vectors pointing to its neighbours add up to a large resultant. A buried
        cell is surrounded, so those unit vectors cancel and the resultant is
        near zero. That resultant length (an isotropy score) is a cheap, robust
        proxy for distance to the front, with no hull or contour needed.

        Growth is ramped from zero in the interior to full at the rim, which is
        what keeps the buried core from churning: once a cell is surrounded it
        stops elongating and dividing, so it generates no new pushing forces.
        """
        cfg = self.cfg
        n = self.n
        if n < 3:
            return np.ones(n), np.zeros(n, dtype=bool)
        pts = np.column_stack((self.x, self.y))
        tree = cKDTree(pts)
        pairs = tree.query_pairs(cfg.front_radius, output_type="ndarray")
        sumx = np.zeros(n)
        sumy = np.zeros(n)
        cnt = np.zeros(n)
        if pairs.size:
            i = pairs[:, 0]
            j = pairs[:, 1]
            v = pts[j] - pts[i]
            d = np.maximum(np.hypot(v[:, 0], v[:, 1]), 1e-9)
            ux = v[:, 0] / d
            uy = v[:, 1] / d
            np.add.at(sumx, i, ux)
            np.add.at(sumy, i, uy)
            np.add.at(cnt, i, 1.0)
            np.add.at(sumx, j, -ux)
            np.add.at(sumy, j, -uy)
            np.add.at(cnt, j, 1.0)
        safe = np.maximum(cnt, 1.0)
        res = np.hypot(sumx, sumy) / safe          # 0 buried .. ~1 lone edge cell
        # cells with very few neighbours are unambiguously at the front
        res = np.where(cnt <= 2, 1.0, res)
        factor = np.clip((res - cfg.front_lo) / (cfg.front_hi - cfg.front_lo),
                         0.0, 1.0)
        # Smooth the growth factor over the neighbourhood. A lone cell that
        # bulges out has fewer neighbours and would otherwise grow faster and
        # bulge more (a fingering instability). Averaging with neighbours damps
        # that runaway, so the rim advances as a smooth, near circular front.
        if pairs.size and cfg.front_smooth_iters > 0:
            i = pairs[:, 0]
            j = pairs[:, 1]
            for _ in range(cfg.front_smooth_iters):
                acc = factor.copy()
                num = np.ones(n)
                np.add.at(acc, i, factor[j])
                np.add.at(acc, j, factor[i])
                np.add.at(num, i, 1.0)
                np.add.at(num, j, 1.0)
                factor = acc / num
        frozen = np.zeros(n, dtype=bool)
        if cfg.freeze_core:
            frozen = (res < cfg.freeze_res) & (cnt >= cfg.freeze_count)
        return factor, frozen

    # -- mechanics ----------------------------------------------------
    def relax(self, iters=None):
        """Resolve overlaps with overdamped position based contact relaxation.

        For each overlapping pair we find the closest points on the two spines,
        push the cells apart along the contact normal, and apply the correction
        at the contact point so off centre contacts also rotate the rods. That
        torque is what makes rods align side by side into nematic domains, which
        point particles never do.
        """
        cfg = self.cfg
        iters = cfg.relax_iters if iters is None else iters
        if self.n < 2:
            return
        reach = cfg.L_div + 2 * cfg.R + 0.5
        mind = 2 * cfg.R
        for _ in range(iters):
            frozen = getattr(self, "_frozen", None)
            if frozen is not None and frozen.shape[0] == self.n and frozen.any():
                x0 = self.x.copy()
                y0 = self.y.copy()
                th0 = self.th.copy()
            pts = np.column_stack((self.x, self.y))
            tree = cKDTree(pts)
            pairs = tree.query_pairs(reach, output_type='ndarray')
            if pairs.size == 0:
                break
            i = pairs[:, 0]
            j = pairs[:, 1]
            p, q = self.endpoints()
            c1, c2 = _segment_closest_batch(p[i], q[i], p[j], q[j])
            d = c2 - c1
            dist = np.hypot(d[:, 0], d[:, 1])
            active = (dist < mind)
            if not active.any():
                continue
            i, j = i[active], j[active]
            c1, c2 = c1[active], c2[active]
            d, dist = d[active], dist[active]
            safe = np.maximum(dist, 1e-9)
            n_hat = d / safe[:, None]
            overlap = mind - dist
            corr = 0.5 * overlap * cfg.relax_stiff

            dx = np.zeros(self.n)
            dy = np.zeros(self.n)
            dth = np.zeros(self.n)
            np.add.at(dx, i, -corr * n_hat[:, 0])
            np.add.at(dy, i, -corr * n_hat[:, 1])
            np.add.at(dx, j, corr * n_hat[:, 0])
            np.add.at(dy, j, corr * n_hat[:, 1])
            # off centre contact -> torque (2D cross product of lever and normal)
            ri = c1 - pts[i]
            rj = c2 - pts[j]
            ti = ri[:, 0] * (-n_hat[:, 1]) - ri[:, 1] * (-n_hat[:, 0])
            tj = rj[:, 0] * (n_hat[:, 1]) - rj[:, 1] * (n_hat[:, 0])
            np.add.at(dth, i, cfg.torque_gain * corr * ti)
            np.add.at(dth, j, cfg.torque_gain * corr * tj)

            self.x = self.x + dx
            self.y = self.y + dy
            self.th = self.th + dth

            # hold the jammed core still: a buried cell acts as a rigid wall, so
            # front pressure does not propagate inward as visible jitter
            frozen = getattr(self, "_frozen", None)
            if frozen is not None and frozen.shape[0] == self.n and frozen.any():
                self.x[frozen] = x0[frozen]
                self.y[frozen] = y0[frozen]
                self.th[frozen] = th0[frozen]

    def compact(self):
        """Drop dead cells so arrays stay small (used after predation)."""
        if self.alive.all():
            return
        m = self.alive
        self.x = self.x[m]
        self.y = self.y[m]
        self.th = self.th[m]
        self.L = self.L[m]
        self.sp = self.sp[m]
        self.lin = self.lin[m]
        self.alive = self.alive[m]
        self._frozen = None
