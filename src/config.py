"""
config.py
=========

One place for every tunable number, so a learner can change a single value and
see the consequence. Defaults are chosen to give a clear range expansion in a
few dozen frames while staying cheap enough to run on a laptop.
"""

from dataclasses import dataclass


@dataclass
class Config:
    # grid and time
    N: int = 200            # diffusion grid is N x N nodes
    dx: float = 1.0         # node spacing (arbitrary length units)
    dt: float = 0.10        # macro time step
    n_frames: int = 40      # snapshots saved per simulation
    steps_per_frame: int = 4

    # seeding (central mixed disk -> range expansion)
    seed_radius: float = 6.0
    n_seed: int = 70
    frac_strain1: float = 0.5

    # rod mechanics
    R: float = 0.55             # rod half width (cap radius), fixed
    L_birth: float = 1.8        # cylinder length at birth
    L_div: float = 3.6          # divide past this length
    div_angle_noise: float = 0.10
    max_cells: int = 1600
    relax_iters: int = 12
    relax_stiff: float = 1.0
    torque_gain: float = 0.30

    # range-expansion front gating
    # Only cells near the growing edge have substrate access and room to grow.
    # Buried cells stop elongating and the packed core freezes, exactly the
    # nutrient-limited range expansion picture: all the action is at the rim.
    front_radius: float = 6.0   # neighbourhood used to judge how buried a cell is
    front_lo: float = 0.08      # isotropy below this -> fully interior (no growth)
    front_hi: float = 0.30      # isotropy above this -> full front growth
    front_smooth_iters: int = 3 # average growth over neighbours: suppress fingering
    freeze_core: bool = True    # hold the jammed core still during contact relax
    freeze_res: float = 0.06    # deep-core isotropy threshold for freezing
    freeze_count: int = 7       # need at least this many neighbours to be frozen

    # growth kinetics
    g_max: float = 2.5          # max elongation (length per time)
    K_S: float = 0.15           # half saturation, primary substrate
    K_M: float = 0.12           # half saturation, secreted metabolite
    K_P: float = 0.14           # half saturation, secreted public good
    Ki: float = 0.16            # toxin concentration for half inhibition
    K_prey: float = 2.0         # prey neighbours for half max predation

    # active vs passive interactions (the paper's central distinction)
    # passive: a by-product leaks for free. active: the cell pays to make the
    # mediating molecule, so the producer's own growth is reduced by a cost.
    cost_public_good: float = 0.30  # growth penalty A pays to make a public good
    pg_gain: float = 1.5            # benefit the public good confers on B
    cost_obligate: float = 0.12     # growth penalty paid in obligate mutualism
    fac_base: float = 0.55          # facultative floor in passive mutualism

    # stoichiometry
    Y_consume: float = 0.060    # substrate or metabolite removed per unit growth
    Y_produce: float = 0.055    # metabolite leaked per unit growth
    Y_toxin: float = 0.160      # toxin produced per unit growth
    active_cost: float = 0.05   # small growth penalty for making a toxin

    # field transport
    D_S: float = 6.0            # primary substrate diffusivity
    D_M: float = 4.0            # metabolite diffusivity
    D_P: float = 4.0            # public good diffusivity
    D_T: float = 5.0            # toxin diffusivity
    S0: float = 1.0             # substrate reservoir at the border
    M_basal: float = 0.07       # small background metabolite to bootstrap mutualism
    decay_M: float = 0.02
    decay_P: float = 0.02
    decay_T: float = 0.03

    # predation
    g_pred: float = 2.0         # predator max elongation from consuming prey
    p_kill: float = 0.022       # per step probability a contacted prey is lysed
    predation_reach: float = 2.4

    seed: int = 7

    @property
    def center(self):
        return self.N * self.dx / 2.0
