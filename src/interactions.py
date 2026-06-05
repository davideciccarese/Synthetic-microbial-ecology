"""
interactions.py
===============

The interaction archetypes of Dolinsek, Goldschmidt and Johnson, "Synthetic
microbial ecology and the dynamic interplay between microbial genotypes", FEMS
Microbiology Reviews 40:961-979 (2016).

That review organises pairwise interactions along two axes that we follow here:

  passive vs active   In a passive interaction the mediating molecule is an
                      inadvertent by-product that leaks for free. In an active
                      interaction the cell pays a metabolic cost to make the
                      molecule, so the producer's own growth is reduced.
  one way vs two way  A unidirectional interaction has one producer and one
                      responder; a bidirectional one is reciprocal.

The review also stresses that competition is pervasive: two genotypes sharing
space and a substrate compete even while one helps the other, so few positive
interactions are purely positive.

Two strains share one space and one set of diffusing chemicals. Strain A (drawn
blue) and strain B (drawn green). Only the chemical coupling changes between
panels; everything spatial is emergent.

Each interaction exposes:
  name, row   a label and a short tag
  signs       effect of the partner on (A, B): +, 0 or -
  kind        "passive" or "active" (or "" where it does not apply)
  seed_frac   strain-B fraction in the founder disk (None = use the global value)
  fields()    the ScalarField set it needs
  step(col, F, cfg, dt) -> mu     per cell elongation rate, after editing fields
"""

import numpy as np
from scipy.spatial import cKDTree

from field import ScalarField


def monod(c, k):
    return c / (k + c)


class Interaction:
    name = "base"
    row = "base"
    signs = ("0", "0")
    kind = ""
    seed_frac = None
    blurb = ""

    def fields(self, cfg):
        return {}

    def step(self, col, F, cfg, dt):
        raise NotImplementedError


def _S(cfg):
    return ScalarField(cfg.N, cfg.dx, cfg.D_S, cfg.dt, c0=cfg.S0,
                       boundary="reservoir", reservoir=cfg.S0)


# ----------------------------------------------------------------------
# 1. Passive unidirectional positive  (commensalism)
#    A grows on substrate S and a by-product M leaks from it for free. B is a
#    high-affinity scavenger of M and occupies a separate niche, so A is
#    unaffected. This is acetate-style cross-feeding (Rosenzweig et al. 1994;
#    Bernstein, Paulson and Carlson 2012).
# ----------------------------------------------------------------------
class CommensalismPassive(Interaction):
    name = "Commensalism (passive, one-way +)"
    row = "commensalism"
    signs = ("0", "+")
    kind = "passive"
    blurb = "A by-product leaks from A for free; B scavenges it. B benefits, A is unaffected."

    def fields(self, cfg):
        return {"S": _S(cfg),
                "M": ScalarField(cfg.N, cfg.dx, cfg.D_M, cfg.dt, c0=0.0,
                                 boundary="zeroflux", decay=cfg.decay_M)}

    def step(self, col, F, cfg, dt):
        sp = col.spine_samples()
        S = F["S"].sample_cells(sp)
        M = F["M"].sample_cells(sp)
        a = col.sp == 0
        b = ~a
        mu = np.zeros(col.n)
        mu[a] = cfg.g_max * monod(S[a], cfg.K_S)
        # high-affinity scavenger of the dilute by-product (small half saturation)
        mu[b] = cfg.g_max * monod(M[b], 0.4 * cfg.K_M)
        dS = np.where(a, -cfg.Y_consume * mu, 0.0) * dt
        dM = np.where(a, 1.8 * cfg.Y_produce * mu, -cfg.Y_consume * mu) * dt
        F["S"].add_cells(sp, dS)
        F["M"].add_cells(sp, dM)
        return mu


# ----------------------------------------------------------------------
# 2. Active unidirectional positive  (a costly public good)
#    Both strains grow on S. A actively secretes a public good P (e.g. an
#    extracellular enzyme such as a beta-lactamase) at a metabolic cost, so A
#    grows more slowly. B exploits P without paying, so B is promoted while A is
#    reduced: the cost of cooperation and the advantage of the non-producer.
# ----------------------------------------------------------------------
class PublicGoodActive(Interaction):
    name = "Public good (active, one-way +)"
    row = "public_good"
    signs = ("-", "+")
    kind = "active"
    blurb = "A pays to secrete a public good that B exploits for free; A is reduced, B is promoted."

    def fields(self, cfg):
        return {"S": _S(cfg),
                "P": ScalarField(cfg.N, cfg.dx, cfg.D_P, cfg.dt, c0=0.0,
                                 boundary="zeroflux", decay=cfg.decay_P)}

    def step(self, col, F, cfg, dt):
        sp = col.spine_samples()
        S = F["S"].sample_cells(sp)
        P = F["P"].sample_cells(sp)
        a = col.sp == 0
        b = ~a
        base = cfg.g_max * monod(S, cfg.K_S)
        mu = np.zeros(col.n)
        mu[a] = base[a] * (1.0 - cfg.cost_public_good)          # A pays the cost
        mu[b] = base[b] * (1.0 + cfg.pg_gain * monod(P[b], cfg.K_P))  # B free-rides
        dS = -cfg.Y_consume * mu * dt
        # only A secretes the public good, in proportion to its (pre-cost) effort
        dP = np.where(a, cfg.Y_produce * base, 0.0) * dt
        F["S"].add_cells(sp, dS)
        F["P"].add_cells(sp, dP)
        return mu


# ----------------------------------------------------------------------
# 3. Passive bidirectional positive  (facultative mutualism)
#    Each strain grows on S on its own, and each leaks a by-product for free that
#    gives the other a growth bonus. Reciprocal but facultative, so the colony
#    still expands while staying intermixed.
# ----------------------------------------------------------------------
class MutualismPassive(Interaction):
    name = "Facultative mutualism (passive +/+)"
    row = "mutualism_passive"
    signs = ("+", "+")
    kind = "passive"
    blurb = "Each strain grows on S alone and gets a free by-product bonus from the other."

    def fields(self, cfg):
        return {"S": _S(cfg),
                "Ma": ScalarField(cfg.N, cfg.dx, cfg.D_M, cfg.dt, c0=0.0,
                                  boundary="zeroflux", decay=cfg.decay_M),
                "Mb": ScalarField(cfg.N, cfg.dx, cfg.D_M, cfg.dt, c0=0.0,
                                  boundary="zeroflux", decay=cfg.decay_M)}

    def step(self, col, F, cfg, dt):
        sp = col.spine_samples()
        S = F["S"].sample_cells(sp)
        Ma = F["Ma"].sample_cells(sp)
        Mb = F["Mb"].sample_cells(sp)
        a = col.sp == 0
        b = ~a
        fb = cfg.fac_base
        s = monod(S, cfg.K_S)
        mu = np.zeros(col.n)
        mu[a] = cfg.g_max * s[a] * (fb + (1 - fb) * monod(Mb[a], cfg.K_M))
        mu[b] = cfg.g_max * s[b] * (fb + (1 - fb) * monod(Ma[b], cfg.K_M))
        dS = -cfg.Y_consume * mu * dt
        dMa = np.where(a, cfg.Y_produce * mu, -0.5 * cfg.Y_consume * mu) * dt
        dMb = np.where(b, cfg.Y_produce * mu, -0.5 * cfg.Y_consume * mu) * dt
        F["S"].add_cells(sp, dS)
        F["Ma"].add_cells(sp, dMa)
        F["Mb"].add_cells(sp, dMb)
        return mu


# ----------------------------------------------------------------------
# 4. Active bidirectional positive  (obligate mutualism)
#    Engineered auxotroph-style cross-feeding. A needs metabolite Mb from B and
#    pays to make Ma; B needs Ma and pays to make Mb. Neither grows without the
#    other close by, so the strains must stay finely intermixed and the colony
#    is the smallest of the positive cases.
# ----------------------------------------------------------------------
class MutualismActive(Interaction):
    name = "Obligate mutualism (active +/+)"
    row = "mutualism_active"
    signs = ("+", "+")
    kind = "active"
    blurb = "Obligate cross-feeding: each strain depends on, and pays to make, the other's metabolite."

    def fields(self, cfg):
        return {"Ma": ScalarField(cfg.N, cfg.dx, cfg.D_M, cfg.dt, c0=cfg.M_basal,
                                  boundary="zeroflux", decay=cfg.decay_M),
                "Mb": ScalarField(cfg.N, cfg.dx, cfg.D_M, cfg.dt, c0=cfg.M_basal,
                                  boundary="zeroflux", decay=cfg.decay_M)}

    def step(self, col, F, cfg, dt):
        sp = col.spine_samples()
        Ma = F["Ma"].sample_cells(sp)
        Mb = F["Mb"].sample_cells(sp)
        a = col.sp == 0
        b = ~a
        c = 1.0 - cfg.cost_obligate
        mu = np.zeros(col.n)
        mu[a] = cfg.g_max * monod(Mb[a], cfg.K_M) * c   # A obligately needs Mb
        mu[b] = cfg.g_max * monod(Ma[b], cfg.K_M) * c   # B obligately needs Ma
        # each makes the metabolite the other needs; each consumes what it needs
        dMa = (np.where(a, cfg.Y_produce * mu, 0.0)
               - np.where(b, cfg.Y_consume * mu, 0.0)) * dt
        dMb = (np.where(b, cfg.Y_produce * mu, 0.0)
               - np.where(a, cfg.Y_consume * mu, 0.0)) * dt
        F["Ma"].add_cells(sp, dMa)
        F["Mb"].add_cells(sp, dMb)
        return mu


# ----------------------------------------------------------------------
# 5. Negative bidirectional  (passive resource competition)
#    Both strains grow on the same S and nothing else. Competition is the
#    inadvertent, passive consequence of a shared resource and demixes the
#    colony into single-strain sectors.
# ----------------------------------------------------------------------
class Competition(Interaction):
    name = "Competition (passive -/-)"
    row = "competition"
    signs = ("-", "-")
    kind = "passive"
    blurb = "Both draw on one shared substrate and nothing else: pure resource competition."

    def fields(self, cfg):
        return {"S": _S(cfg)}

    def step(self, col, F, cfg, dt):
        sp = col.spine_samples()
        S = F["S"].sample_cells(sp)
        mu = cfg.g_max * monod(S, cfg.K_S)
        F["S"].add_cells(sp, -cfg.Y_consume * mu * dt)
        return mu


# ----------------------------------------------------------------------
# 6. Active negative unidirectional  (amensalism via an antibiotic)
#    A actively makes a diffusing antibiotic T at a small cost. T inhibits B; B
#    has essentially no effect back on A, so B is excluded from a halo around the
#    A patches.
# ----------------------------------------------------------------------
class Amensalism(Interaction):
    name = "Amensalism (active, one-way -)"
    row = "amensalism"
    signs = ("0", "-")
    kind = "active"
    blurb = "A makes a diffusing antibiotic at a small cost that inhibits B; B has no effect back."

    def fields(self, cfg):
        return {"S": _S(cfg),
                "T": ScalarField(cfg.N, cfg.dx, cfg.D_T, cfg.dt, c0=0.0,
                                 boundary="zeroflux", decay=cfg.decay_T)}

    def step(self, col, F, cfg, dt):
        sp = col.spine_samples()
        S = F["S"].sample_cells(sp)
        T = F["T"].sample_cells(sp)
        a = col.sp == 0
        b = ~a
        base = cfg.g_max * monod(S, cfg.K_S)
        mu = np.zeros(col.n)
        mu[a] = base[a] * (1.0 - cfg.active_cost)            # small cost to make T
        mu[b] = base[b] * (1.0 / (1.0 + (T[b] / cfg.Ki) ** 2))
        dS = -cfg.Y_consume * mu * dt
        dT = np.where(a, cfg.Y_toxin * base, 0.0) * dt
        F["S"].add_cells(sp, dS)
        F["T"].add_cells(sp, dT)
        kill_b = b & (T > 3.0 * cfg.Ki)
        if kill_b.any():
            roll = col.rng.random(col.n)
            col.alive[kill_b & (roll < 0.10)] = False
        return mu


# ----------------------------------------------------------------------
# 7. Predation / parasitism
#    Prey A grows on S. Predator B grows only by contacting prey, and a contacted
#    prey cell is lysed at a low per step rate. Starting from few predators, the
#    predator advances as green inroads into a growing prey colony.
# ----------------------------------------------------------------------
class Predation(Interaction):
    name = "Predation / parasitism (-/+)"
    row = "predation"
    signs = ("-", "+")
    kind = ""
    seed_frac = 0.18     # start from a prey colony with a minority of predators

    blurb = "Predator B grows on contact with prey A and lyses it; it tracks the prey it consumes."

    def fields(self, cfg):
        return {"S": _S(cfg)}

    def step(self, col, F, cfg, dt):
        sp = col.spine_samples()
        S = F["S"].sample_cells(sp)
        a = col.sp == 0           # prey
        b = ~a                    # predator
        mu = np.zeros(col.n)
        mu[a] = cfg.g_max * monod(S[a], cfg.K_S)         # prey grows on S
        F["S"].add_cells(sp, np.where(a, -cfg.Y_consume * mu, 0.0) * dt)

        pts = np.column_stack((col.x, col.y))
        ia = np.where(a)[0]
        ib = np.where(b)[0]
        if ia.size and ib.size:
            tree = cKDTree(pts[ia])
            counts = np.zeros(ib.size)
            hit_load = np.zeros(ia.size)
            for k, j in enumerate(ib):
                nb = tree.query_ball_point(pts[j], cfg.predation_reach)
                counts[k] = len(nb)
                for m in nb:
                    hit_load[m] += 1.0
            mu[ib] = cfg.g_pred * monod(counts, cfg.K_prey)   # grows from prey
            contacted = ia[hit_load > 0]
            if contacted.size:
                roll = col.rng.random(contacted.size)
                col.alive[contacted[roll < cfg.p_kill]] = False
        return mu


# ordered list, matching the review's progression: positive (passive then
# active, one way then two way), then negative (competition, amensalism,
# predation)
ALL = [
    CommensalismPassive(),
    PublicGoodActive(),
    MutualismPassive(),
    MutualismActive(),
    Competition(),
    Amensalism(),
    Predation(),
]
