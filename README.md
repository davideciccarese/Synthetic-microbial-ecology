# Microbial interaction archetypes in a range-expanding rod colony

A small, readable, spatially explicit model that plays out every binary interaction type in

> Dolinšek J, Goldschmidt F, Johnson DR. Synthetic microbial ecology and the dynamic interplay between microbial genotypes. *FEMS Microbiology Reviews* 40(6):961-979 (2016). https://doi.org/10.1093/femsre/fuw024

The same colony of rod-shaped cells, seeded the same way every time, is run under seven different chemical couplings. Only the coupling changes between panels. Everything spatial (who reaches the front, who ends up next to whom, whether the strains mix or demix) is emergent from local growth and diffusion, not scripted.

![interaction archetypes, coloured by strain](figures/interactions_strain.gif)

The figure above reproduces the logic of Figure 1 of the review as a living colony rather than a set of abundance curves.

## What the model is

Each cell is a spherocylinder (a rod with hemispherical caps) that:

- elongates at a rate set by the local concentration of whatever chemical it needs,
- divides along its long axis once it grows past a threshold length,
- pushes and rotates its neighbours through overdamped contact mechanics, which is what makes the rods align side by side into nematic domains,
- reads and writes a shared set of diffusing chemical fields.

The cells start as a small mixed disk of founders at the centre and grow outward, so the colony is a range expansion. Substrate is supplied from the border, as on an agar plate. Secreted metabolites and antibiotics use zero-flux boundaries and stay local, so they accumulate next to the cells that make them. That is the localising effect of spatial structure the review stresses: in a structured habitat, diffusion keeps secreted molecules concentrated where they are produced, which is exactly what lets cross-feeding and toxin halos matter.

Growth is restricted to the colony rim. A cell judges how buried it is from whether its neighbours surround it, and once it is enclosed it stops elongating and the packed core freezes. This is the nutrient-limited range-expansion picture, where only the outer shell has resource access, and it is also what keeps the founder mosaic in the core fixed while sectors widen at the front. The growth factor is smoothed across neighbours so a single protruding cell cannot run away, which keeps the expanding front close to circular.

The rod mechanics reuse the engine conventions from earlier work on rod-shaped cross-feeding colonies (Ciccarese et al. 2022, see references), and the nematic alignment of growing rods follows Volfson et al. (2008).

## The seven interactions

The review organises pairwise interactions along two axes, which this repository follows directly. The first is **passive versus active**: in a passive interaction the mediating molecule is an inadvertent by-product that leaks for free, whereas in an active interaction the cell pays a metabolic cost to make the molecule, so the producer's own growth is reduced. The second is **unidirectional versus bidirectional** (one-way versus reciprocal). The review also stresses that competition is pervasive, so two genotypes sharing space and a substrate compete even while one helps the other.

Strain A is drawn blue, strain B green. The sign pair is the effect of the partner on each strain.

| Panel | Class | Sign (A, B) | Chemical coupling | What to look for |
|---|---|---|---|---|
| Commensalism | passive, one-way + | (0, +) | A leaks a by-product M for free; B is a high-affinity scavenger of M and does not touch substrate S | A is unaffected and dominates; B persists threaded through it |
| Public good | active, one-way + | (-, +) | A pays a cost to secrete a diffusing public good P; B exploits P without paying | B (the non-producer) overtakes, the classic cost of cooperation |
| Facultative mutualism | passive, two-way + | (+, +) | Each grows on S alone and gets a free by-product bonus from the other | Balanced, intermixed, the colony still expands |
| Obligate mutualism | active, two-way + | (+, +) | Each strain requires, and pays to make, the metabolite the other needs | Tightest intermixing and the smallest colony, the cost of obligate exchange |
| Competition | passive, two-way - | (-, -) | One shared substrate, nothing else | Clean demixing into single-strain sectors at the front |
| Amensalism | active, one-way - | (0, -) | A makes a diffusing antibiotic at a small cost; B is inhibited, B has no effect back | A gains territory, B suppressed at the contact margin |
| Predation / parasitism | exploitation | (-, +) | Predator B grows only on contact with prey A and lyses it | A minority of predators forms green inroads into a growing prey colony |

The four positive panels are deliberately distinct: commensalism is producer-dominated (A unaffected), the public good flips to non-producer-dominated (A pays), facultative mutualism is balanced and expands, and obligate mutualism is balanced but small and finely mixed. This is the active/passive contrast the review introduces.

A second rendering colours each cell by its founding lineage, which makes the range-expansion sectors visible as stripes:

![interaction archetypes, coloured by founder lineage](figures/interactions_lineage.gif)


Note, in line with the review, that no interaction here is purely positive. Positive couplings carry a competitive element at the same time, because the strains still share space and, in most panels, substrate.

## Model and equations

The model couples a few diffusing chemical fields to an individual-based colony of rods. All symbols are tunable in `config.py`.

**Diffusion fields.** Each chemical field $c(\mathbf{x},t)$ obeys a reaction-diffusion equation on a regular grid,

$$\frac{\partial c}{\partial t} = D\,\nabla^2 c \;-\; \lambda\,c \;+\; \sum_{i}\sigma_i\,\mathbb{1}_{\mathbf{x}\in\,\text{cell }i},$$

with diffusivity $D$, first-order decay $\lambda$, and a per-cell source or sink $\sigma_i$ deposited on the grid nodes under each cell body. Boundary conditions encode the habitat: the substrate uses a Dirichlet (reservoir) condition $c=S_0$ on the border, like an agar plate fed from outside, while secreted molecules use a Neumann zero-flux condition $\partial c/\partial n=0$ so they stay local. The field is integrated with explicit Euler and stability sub-stepping such that $D\,\Delta t_\text{sub}/\Delta x^2 < 1/4$.

**Growth kinetics.** Writing the Monod function $m(c,K)=\dfrac{c}{K+c}$, the per-cell elongation rate $\mu_i$ for each interaction is

$$
\begin{aligned}
\text{Commensalism (passive):}\quad & \mu_A = g\,m(S,K_S), \quad \mu_B = g\,m(M,\,0.4K_M)\\
\text{Public good (active):}\quad & \mu_A = g\,m(S,K_S)\,(1-c_\text{pg}), \quad \mu_B = g\,m(S,K_S)\,[\,1+\gamma\,m(P,K_P)\,]\\
\text{Facultative mutualism:}\quad & \mu_A = g\,m(S,K_S)\,[\,f+(1-f)\,m(M_b,K_M)\,], \quad \text{(B symmetric)}\\
\text{Obligate mutualism (active):}\quad & \mu_A = g\,m(M_b,K_M)\,(1-c_\text{ob}), \quad \mu_B = g\,m(M_a,K_M)\,(1-c_\text{ob})\\
\text{Competition:}\quad & \mu_A = \mu_B = g\,m(S,K_S)\\
\text{Amensalism (active):}\quad & \mu_A = g\,m(S,K_S)\,(1-c_T), \quad \mu_B = g\,m(S,K_S)\,\frac{1}{1+(T/K_i)^2}\\
\text{Predation:}\quad & \mu_A = g\,m(S,K_S), \quad \mu_B = g_p\,m(n_\text{prey},K_\text{prey})
\end{aligned}
$$

Active interactions carry an explicit cost ($c_\text{pg}, c_\text{ob}, c_T$) that the producer pays, which is the passive-versus-active distinction the review emphasises. In predation, $n_\text{prey}$ is the number of prey spines within a reach radius of a predator, and each contacted prey cell is lysed with probability $p_\text{kill}$ per step. Each growing cell edits the fields under it: consumption $-Y_c\,\mu$, by-product or public good $+Y_p\,\mu$, antibiotic $+Y_T\,\mu_\text{base}$.

**Range-expansion front.** Only the rim grows. For a cell $i$ with neighbours $j$ within radius $r_f$, an isotropy score

$$\rho_i = \frac{1}{n_i}\left\lVert\,\sum_{j}\hat{\mathbf{u}}_{ij}\right\rVert,\qquad \hat{\mathbf{u}}_{ij}=\frac{\mathbf{x}_j-\mathbf{x}_i}{\lVert\mathbf{x}_j-\mathbf{x}_i\rVert},$$

is near $0$ for a surrounded (buried) cell, whose neighbour directions cancel, and near $1$ for an exposed edge cell. The growth factor $\phi_i=\mathrm{clip}\!\left(\dfrac{\rho_i-\rho_\text{lo}}{\rho_\text{hi}-\rho_\text{lo}},\,0,\,1\right)$, smoothed over neighbours to suppress single-cell fingering, scales elongation as $\dot L_i=\mu_i\,\phi_i$. Deeply buried cells ($\rho_i<\rho_\text{freeze}$ with many neighbours) are frozen in place, which is what keeps the founder mosaic fixed in the core while sectors widen at the front.

**Cell mechanics.** Each cell is a spherocylinder with centre $\mathbf{x}_i$, orientation $\theta_i$, cylinder length $L_i$ and fixed radius $R$, with spine endpoints $\mathbf{x}_i\pm\frac{L_i}{2}(\cos\theta_i,\sin\theta_i)$. A cell divides when $L_i\ge L_\text{div}$ into two daughters of length $\tfrac{1}{2}(L_i-2R)$. Overlaps are resolved by overdamped, position-based contact relaxation: for a contacting pair the closest points on the two spines give a contact normal $\hat{\mathbf{n}}$ and overlap $\delta=2R-d$; each cell is displaced by $\pm\tfrac{1}{2}k\,\delta\,\hat{\mathbf{n}}$ and rotated by a torque proportional to the off-centre lever $\mathbf{r}\times\hat{\mathbf{n}}$, iterated each step. That contact torque is what aligns growing rods into nematic domains (Volfson et al. 2008).

## Install

```bash
git clone <your-fork-url>
cd ecomodel
pip install -r requirements.txt
```

Dependencies are numpy, scipy, matplotlib and pillow. Python 3.9 or newer.

## Run

Build the three multipanel GIFs (the hero figures), cells drawn as capsule contours:

```bash
python scripts/make_multipanel.py --mode strain
python scripts/make_multipanel.py --mode lineage
python scripts/make_fields_multipanel.py
```

The first colours cells by strain, the second by founder lineage (with a colour bar), and the third shows the representative diffusion field each colony shapes as a heatmap with the cells outlined on top.

Run and render a single interaction on its own:

```bash
python scripts/run_single.py competition
python scripts/run_single.py mutualism --mode lineage
python scripts/run_single.py predation
```

Render the diffusion fields a colony shapes, the substrate it depletes and the metabolites or antibiotic it builds up locally:

```bash
python scripts/make_fields.py --row commensalism
python scripts/make_fields.py --row amensalism
python scripts/make_fields.py --row mutualism
```

![cells next to the fields they shape](figures/fields_commensalism.gif)

The left panel is the colony, coloured by strain. Each further panel is one chemical field as a heatmap with the colony drawn on top. Substrate is supplied from the border and gets drawn down under the colony, while secreted molecules stay concentrated next to the cells that make them. That localisation is the spatial-structure effect the review emphasises.

Valid interaction tags: `commensalism`, `commensalism_compS`, `commensalism_compSM`, `mutualism`, `competition`, `predation`, `amensalism`.

Output GIFs land in `figures/`. A full multipanel run takes well under a minute on a laptop.

## Repository layout

```
ecomodel/
  src/
    config.py        every tunable number in one dataclass
    rods.py          spherocylinder mechanics: elongation, division, contact relaxation, lineage
    field.py         scalar diffusion fields with reservoir or zero-flux boundaries
    interactions.py  the seven archetypes, each as a small growth-plus-exchange rule
    sim.py           wires colony, fields and one interaction, captures snapshots
    render.py        strain and lineage colouring, panel drawing, GIF assembly
  scripts/
    make_multipanel.py
    make_fields_multipanel.py
    make_fields.py
    run_single.py
  docs/
    interactions.md  longer notes on each interaction and its link to the review
  figures/           generated GIFs and preview stills
```

## How to use this for teaching

The code is deliberately the minimal shared core of published off-lattice colony engines, so it is a good thing to read before reaching for CellModeller, iDynoMiCS or BSim. Suggested exercises:

1. Open `config.py` and lower the substrate diffusivity `D_S` or raise the consumption yield `Y_consume`. Watch the smooth front go unstable and finger.
2. In `interactions.py`, give one strain a higher `g_max` under competition and watch one sector sweep the whole front.
3. In the obligate mutualism, raise `cost_obligate` or lower `M_basal` and watch the colony shrink and struggle to bootstrap, the honest price of obligate exchange.
4. Raise the predation kill probability `p_kill` and find the point where the predator drives the prey, and itself, to collapse.
5. Add an eighth interaction class. A good target is the beyond-binary rock-paper-scissors antibiotic system from the review (Kerr et al. 2002, Kelsic et al. 2015), which needs three strains and three toxin fields.

Each interaction is a single class with two methods (`fields` and `step`), so adding one is a short, self-contained edit.

## Modelling notes and honest limitations

- The cells are immobile once placed and grow only by division and mechanical pushing, which is the right picture for cells held in a gel or packed in a colony. Biomass never diffuses.
- Growth laws are Monod, with explicit costs for the active interactions and a high-affinity term for the commensal scavenger. Stoichiometry is a single yield per unit elongation, not a full metabolic network.
- Predation is contact-based (Bdellovibrio-like), not diffusible.
- The numbers are tuned for clear teaching images in a few dozen frames, not fitted to a specific organism. They are a starting point for your own parameterisation.

## References

The interaction framework, and the source of the seven archetypes:

> Dolinšek J, Goldschmidt F, Johnson DR. Synthetic microbial ecology and the dynamic interplay between microbial genotypes. *FEMS Microbiology Reviews* 40(6):961-979 (2016). https://doi.org/10.1093/femsre/fuw024

Examples and concepts the panels draw on, as cited within the review:

> Rosenzweig RF, Sharp RR, Treves DS, Adams J. Microbial evolution in a simple unstructured environment: genetic differentiation in *Escherichia coli*. *Genetics* 137:903-917 (1994).

> Bernstein HC, Paulson SD, Carlson RP. Synthetic *Escherichia coli* consortia engineered for syntrophy demonstrate enhanced biomass productivity. *Journal of Biotechnology* 157:159-166 (2012). https://doi.org/10.1016/j.jbiotec.2011.10.001

> Morris JJ, Lenski RE, Zinser ER. The Black Queen Hypothesis: evolution of dependencies through adaptive gene loss. *mBio* 3:e00036-12 (2012). https://doi.org/10.1128/mBio.00036-12

> Hallatschek O, Hersen P, Ramanathan S, Nelson DR. Genetic drift at expanding frontiers promotes gene segregation. *PNAS* 104:19926-19930 (2007). https://doi.org/10.1073/pnas.0710150104

> Müller MJI, Neugeboren BI, Nelson DR, Murray AW. Genetic drift opposes mutualism during spatial population expansion. *PNAS* 111(3):1037-1042 (2014). https://doi.org/10.1073/pnas.1313285111

Growth kinetics and rod mechanics:

> Monod J. The growth of bacterial cultures. *Annual Review of Microbiology* 3:371-394 (1949). https://doi.org/10.1146/annurev.mi.03.100149.002103

> Volfson D, Cookson S, Hasty J, Tsimring LS. Biomechanical ordering of dense cell populations. *PNAS* 105(40):15346-15351 (2008). https://doi.org/10.1073/pnas.0706805105

> Rudge TJ, Steiner PJ, Phillips A, Haseloff J. Computational modeling of synthetic microbial biofilms. *ACS Synthetic Biology* 1(8):345-352 (2012). https://doi.org/10.1021/sb300031n

The rod-mechanics engine conventions reuse those validated in:

> Ciccarese D, Micali G, Borer B, Ruan C, Or D, Johnson DR. Rare and localized events stabilize microbial community composition and patterns of spatial self-organization in a fluctuating environment. *The ISME Journal* 16:1453-1463 (2022). https://doi.org/10.1038/s41396-022-01189-9

The closest-points-between-segments routine follows Ericson C, *Real-Time Collision Detection*, Morgan Kaufmann (2005).

## License

MIT. See [LICENSE](LICENSE).
